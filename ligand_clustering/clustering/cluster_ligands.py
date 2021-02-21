# coding: utf-8

import os
import ast
import shutil
import logging
import pandas as pd
import itertools

from rdkit import Chem
from rdkit.Chem import Draw
from utils import DbScript
from populate_ligands import PopulateLigandsScript


class ClusteringScript(DbScript):
    CONNECTIVITY_MATCHES_DB = "matching_connectivity"
    CONNECTIVITY_CLUSTERS_DB = "clusters_connectivity"
    CONNECTIVITY_MAPPING_DB = "mapping_connectivity"
    GEOMETRY_MATCHES_DB = "matching_geometry_and_weight"
    GEOMETRY_CLUSTERS_DB = "clusters_geometry_and_weight"
    GEOMETRY_MAPPING_DB = "mapping_geometry_and_weight"

    def __init__(self, draw_ligands=True, logging_level=logging.INFO):
        super().__init__("Cluster Connectivity", logging_level=logging_level,
                         tables_to_drop=[self.CONNECTIVITY_MATCHES_DB, self.CONNECTIVITY_CLUSTERS_DB,
                                         self.CONNECTIVITY_MAPPING_DB, self.GEOMETRY_MATCHES_DB,
                                         self.GEOMETRY_CLUSTERS_DB, self.GEOMETRY_MAPPING_DB])
        self.bond_pattern_sets = [[Chem.MolFromSmarts('[!D1]#[!D1]')],
                                  [Chem.MolFromSmarts('[!D1]=*=[!D1]')],
                                  [Chem.MolFromSmarts('[D1]#[D2]-*'), Chem.MolFromSmarts('[D1]=[D2]=*')]]
        self.draw_ligands = draw_ligands
        self.ligands_df = None
        self.counts_df = None
        self.ligands_dict = None
        self.pattern_dict = None

    def run(self):
        self.logger.info("Running clustering script")

        super().run()
        self.counts_df = self._get_counts()
        self.ligands_df = self._get_ligands()
        self.ligands_dict = self.ligands_df.set_index("pdbid")['molecule'].to_dict()
        self.pattern_dict = self.ligands_df.set_index("pdbid")['pattern'].to_dict()

        connectivity_df = self.match_connectivity()
        geometry_df = self.match_geometry_and_weight(connectivity_df)
        if self.draw_ligands:
            self.draw_ligand_clusters(self.GEOMETRY_CLUSTERS_DB)

        self.logger.info("Finished clustering script")

    def match_connectivity(self):
        self.logger.info("Matching ligands by connectivity")

        matches = [(pdbid, pdbid) for pdbid in self.ligands_df.pdbid]

        for name, group_df in \
                self.ligands_df.loc[self.ligands_df.atoms > 0, :].groupby(["atoms", "rings"]):
            self.logger.debug("Group %s (%s)", name, group_df.shape[0])
            group_matches = 0
            number_of_atoms = name[0]
            pair_count = self._get_number_of_pairs(group_df.shape[0])

            for pdbid1, pdbid2 in itertools.combinations(group_df.pdbid, 2):
                self.logger.debug("(%s, %s)", pdbid1, pdbid2)
                if number_of_atoms == 1 or \
                        self.ligands_dict[pdbid1].HasSubstructMatch(self.pattern_dict[pdbid2], useChirality=False):
                    matches.append((pdbid1, pdbid2))
                    matches.append((pdbid2, pdbid1))
                    group_matches += 1

            self.logger.debug("%s out of %s pairs matched", group_matches, pair_count)

        connectivity_df = self.format_and_cluster(matches, self.CONNECTIVITY_MATCHES_DB, self.CONNECTIVITY_CLUSTERS_DB, self.CONNECTIVITY_MAPPING_DB)

        return connectivity_df

    def match_geometry_and_weight(self, prior_matches_df):
        self.logger.info("Matching ligands by geometry")
        matches = [(pdbid, pdbid) for pdbid in prior_matches_df.index.unique()]
        checked = set()

        for pdbid1, pdbid2 in prior_matches_df.itertuples():
            pair = tuple(sorted((pdbid1, pdbid2)))
            atom_num = Chem.rdchem.Mol.GetNumHeavyAtoms(self.ligands_dict[pdbid1])

            if pair in checked or pdbid1 == pdbid2:
                continue
            elif atom_num < 3:
                is_match = self._check_weight_substitutes(self.ligands_dict[pdbid1], self.ligands_dict[pdbid2])
            else:
                m1 = self.ligands_dict[pdbid1]
                m2 = self.ligands_dict[pdbid2]
                m2_smart = self.pattern_dict[pdbid2]
                potential_matches = m1.GetSubstructMatches(m2_smart, uniquify=False, useChirality=False)
                if potential_matches:
                    self._fix_chirality_tags([m1, m2])

                for match in potential_matches:
                    try:
                        m1_renumbered = Chem.RenumberAtoms(m1, match)
                        m1_renumbered.UpdatePropertyCache()
                        is_match = self._is_chirality_match(m1_renumbered, m2)
                    except:
                        print(pdbid1, pdbid2)
                        is_match = False

                    if is_match:
                        is_match = self._check_bonds(m1_renumbered, m2) and \
                                   self._check_weight_substitutes(m1_renumbered, m2)
                        break

            if is_match:
                matches.append((pdbid1, pdbid2))
                if pdbid1 != pdbid2: matches.append((pdbid2, pdbid1))

            checked.add(pair)

        geometry_df = self.format_and_cluster(matches, self.GEOMETRY_MATCHES_DB, self.GEOMETRY_CLUSTERS_DB, self.GEOMETRY_MAPPING_DB)
        return geometry_df

    def format_and_cluster(self, matches, matches_table, clusters_table, mapping_table, clustering_must_be_crisp=True):
        matches_df = pd.DataFrame(matches, columns=["pdbid", "match"])
        matches_df = matches_df.set_index('pdbid')
        matches_df = matches_df.sort_index()

        self.save_result_df(matches_df, matches_table)
        self.cluster_ligand_matches(matches_df, clusters_table, mapping_table, clustering_must_be_crisp)

        return matches_df

    def cluster_ligand_matches(self, matching_df, output_name, mapping_table, clustering_must_be_crisp):
        self.logger.info("Clustering data into {0}".format(output_name))
        matching_df = matching_df.join(self.counts_df, how="left")
        mapping_df = pd.DataFrame()

        grouped_df = matching_df.groupby(["pdbid"], as_index=False).agg(
            cluster=("match", frozenset), instances=("instances", "first"), deposits=("deposits", "first"))
        grouped_df = grouped_df.groupby(["cluster"]).agg(
            instances=("instances", "sum"), deposits=("deposits", "sum"), length=("instances", "count")).reset_index()
        grouped_df.loc[:, "cluster"] = grouped_df.cluster.apply(self._get_cluster_ligand_order)
        grouped_df.loc[:, "name"] = grouped_df.cluster.apply(self._get_cluster_name)
        grouped_df = grouped_df.sort_values("instances", ascending=False).reset_index(drop=True)

        mapping_dfs = []
        grouped_df.apply(self._create_mapping, axis=1, mapping_dfs=mapping_dfs)
        mapping_df = pd.concat(mapping_dfs, ignore_index=True)

        if clustering_must_be_crisp:
            self._perform_sanity_check(grouped_df)

        self.save_result_df(grouped_df, output_name)
        self.save_result_df(mapping_df, mapping_table)

        return grouped_df

    def draw_ligand_clusters(self, cluster_db_name, only_multiple=True, output_folder="data", useSVG=False):
        self.logger.info("Drawing ligands from %s", cluster_db_name)
        cluster_df = pd.read_sql_query("select * from {table};".format(table=cluster_db_name), self.conn)

        output_folder = os.path.join(output_folder, "images_for_" + cluster_db_name)
        if os.path.exists(output_folder):
            shutil.rmtree(output_folder)
        os.makedirs(output_folder)

        for row in cluster_df.itertuples():
            ligand_names = [n.strip() for n in row.cluster.strip('{}').split(',')]

            if not only_multiple or len(ligand_names) > 1:
                ligands = [self.ligands_dict[pdbid] for pdbid in ligand_names]
                image_filename = os.path.join(output_folder, "cluster_{0}_[{1}].png".format(row.index, "_".join(ligand_names)))

                try:
                    with Draw.MolsToGridImage(ligands, legends=ligand_names, molsPerRow=min(4, len(ligands)),
                                              useSVG=useSVG) as img:
                        img.save(image_filename, "PNG")
                except:
                    self.logger.warning("Could not create image for cluster %s [%s]", row.index, "_".join(ligand_names))

    def _perform_sanity_check(self, grouped_df):
        ligand_list = grouped_df.cluster.apply(pd.Series).stack().reset_index(drop=True)
        duplicates = ligand_list.loc[ligand_list.duplicated()].unique()
        if len(duplicates) > 0:
            self.logger.error("Ligands in mutliple clusters: %s", duplicates)
            for pdbid in duplicates:
                self.logger.info("  %s: SMILES: %s, SMARTS: %s",
                                 pdbid,
                                 Chem.MolToSmiles(self.ligands_dict[pdbid]),
                                 Chem.MolToSmarts(self.pattern_dict[pdbid])
                )

    def _fix_chirality_tags(self, molecules):
        for molecule in molecules:
            self._fix_PS_atom_chirality_tag(molecule)

    def _fix_PS_atom_chirality_tag(self, m1):
        fragments_to_fix = m1.GetSubstructMatches(Chem.MolFromSmarts("[P,S](O)(*)(*)=O"), uniquify=False,
                                                  useChirality=False)
        for fragment in fragments_to_fix:
            m1.GetAtoms()[fragment[0]].SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

    def _is_chirality_match(self, m1, m2):
        for atom1, atom2 in zip(m1.GetAtoms(), m2.GetAtoms()):
            tag1 = atom1.GetChiralTag()
            tag2 = atom2.GetChiralTag()
            self.logger.debug("(1) %s, %s (2) %s, %s", atom1.GetIdx(), tag1, atom2.GetIdx(), tag2)

            if tag1 == Chem.ChiralType.CHI_OTHER or tag2 == Chem.ChiralType.CHI_OTHER:
                self.logger.error("Unknown chirality tag")

            if tag1 != Chem.ChiralType.CHI_UNSPECIFIED and tag2 != Chem.ChiralType.CHI_UNSPECIFIED:
                self.logger.debug("Both chiral")
                bond_order1 = self._get_bond_order(atom1, tag1)
                bond_order2 = self._get_bond_order(atom2, tag2)
                match_chiral = self._match_rotation(bond_order1, bond_order2)
            elif tag1 == Chem.ChiralType.CHI_UNSPECIFIED and tag2 == Chem.ChiralType.CHI_UNSPECIFIED:
                self.logger.debug("Both unspecified")
                match_chiral = True
            else:
                self.logger.debug("One chiral, one unspecified")
                match_chiral = False

            if not match_chiral:
                return False

        return True

    def _match_rotation(self, bond_order1, bond_order2):
        match_rotation = False

        # Rotate the element of bonds, until match or mismatch
        for i in range(len(bond_order1)):
            bond_order2 = bond_order2[-1:] + bond_order2[:-1]
            self.logger.debug('BO: %s %s', bond_order1, bond_order2)
            if bond_order1 == bond_order2:
                match_rotation = True
                break

        return match_rotation

    def _check_bonds(self, molecule1, molecule2):
        """
        The molecules will be checked against a set of patterns
        If both molecules match this pattern it will be checked
        if the pattern matches in the same place
        If multiple patterns are specified for a set
        then they are treated as equivalent ie. the same place
        in the molecule must match any of the patterns.
        """
        rejected = False

        for patts in self.bond_pattern_sets:
            matches = [(molecule1.GetSubstructMatch(patt, useChirality=False),
                        molecule2.GetSubstructMatch(patt, useChirality=False), patt) for patt in patts]
            empty = True
            pattern_match = False

            for product in itertools.product(matches, repeat=2):
                match1 = product[0][0]
                match2 = product[1][1]
                patt1 = product[0][2]
                patt2 = product[1][2]
                self.logger.debug("%s, %s, %s, %s", patt1, patt2, match1, match2)

                self.logger.debug("%s %s", match1, match2)
                if match1 or match2:
                    empty = False
                if match1 and ((match1 == match2) or (match1 == tuple(reversed(match2)))):
                    pattern_match = True
                    self.logger.debug("Accepted")
                    break

            if empty:
                self.logger.debug("Accepted")
                pattern_match = True

            if not pattern_match:
                self.logger.debug("Rejected")
                rejected = True
                break

        matching = not rejected
        return matching

    def _get_bond_order(self, atoms, tag):
        bond_order1 = []
        for b in atoms.GetBonds():
            bond_order1.append(b.GetOtherAtomIdx(atoms.GetIdx()))
        if tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            bond_order1 = list(reversed(bond_order1))
        return bond_order1

    def _get_ligands(self):
        ligands_df = pd.read_sql_query("""
                select pdbid, mol_send(molecule) as molecule, atoms, rings, aromatic_rings, weight
                from {ligands}
                where molecule is not null
        """.format(ligands=PopulateLigandsScript.LIGANDS_DB), self.conn)

        params = Chem.AdjustQueryParameters()
        params.makeAtomsGeneric = True
        params.makeBondsGeneric = True
        params.adjustRingCount = True

        ligands_df.loc[:, "molecule"] = ligands_df.loc[:, "molecule"].apply(lambda m: Chem.Mol(m.tobytes()))
        ligands_df.loc[:, "pattern"] = ligands_df.loc[:, "molecule"].apply(lambda m: Chem.AdjustQueryProperties(m, params))

        return ligands_df

    def _get_counts(self):
        counts_df =  pd.read_sql_query("""
            select * from {counts}
        """.format(counts=PopulateLigandsScript.COUNTS_DB), self.conn).set_index('pdbid')

        return counts_df

    @staticmethod
    def _get_number_of_pairs(count):
        return int(count*(count - 1)/2)

    def _check_weight_substitutes(self, m1, m2):
        for atom1, atom2 in zip(m1.GetAtoms(), m2.GetAtoms()):
            symbol1 = atom1.GetSymbol()
            symbol2 = atom2.GetSymbol()

            if not (symbol1 == symbol2 or self._get_weight_group(atom1) == self._get_weight_group(atom2)):
                return False

        return True

    @staticmethod
    def _get_weight_group(atom):
        atomic_num = atom.GetAtomicNum()

        if atomic_num < 5:
            return 0
        elif atomic_num < 10:
            return 1
        elif atomic_num < 11:
            return 2
        elif atomic_num < 14:
            return 3
        elif atomic_num < 18:
            return 4
        elif atomic_num < 19:
            return 5
        elif atomic_num < 21:
            return 6
        elif atomic_num < 32:
            return 7
        elif atomic_num < 36:
            return 8
        elif atomic_num < 37:
            return 9
        elif atomic_num < 39:
            return 10
        elif atomic_num < 50:
            return 11
        elif atomic_num < 54:
            return 12
        elif atomic_num < 55:
            return 13
        else:
            return 14

    @staticmethod
    def _get_cluster_name(list_of_ligands):
        if len(list_of_ligands) == 1:
            return list_of_ligands[0]
        else:
            return list_of_ligands[0] + "-like"

    def _get_cluster_ligand_order(self, set_of_ligands):
        ligands = sorted(list(set_of_ligands))

        if len(ligands) == 1:
            return ligands
        else:
            ligand_counts = self.counts_df.loc[self.counts_df.index.isin(ligands), :].instances.sort_index().to_list()
            return [x for _, x in sorted(zip(ligand_counts, ligands), reverse=True)]

    def _create_mapping(self, row, mapping_dfs):
        mapping_dfs.append(pd.DataFrame({"ligand": row["cluster"], "cluster_name": [row["name"]] * len(row.cluster)}))


# Main script #
if __name__ == '__main__':
    script = ClusteringScript(logging_level=logging.INFO)
    script.run()
