# coding: utf-8

import logging
import os
import pandas as pd
import re
import sys

from utils import DbScript
from rdkit import Chem
from rdkit import RDLogger

from ligand_fixes import fixes as LIGAND_FIXES

class PopulateLigandsScript(DbScript):
    RAW_DATA_DB = "raw_data"
    COUNTS_DB = "counts"
    LIGANDS_DB = "ligands"

    def __init__(self, logging_level=logging.INFO):
        super().__init__("Populate", logging_level=logging_level,
                         tables_to_drop=[self.RAW_DATA_DB, self.COUNTS_DB, self.LIGANDS_DB])
        lg = RDLogger.logger()
        lg.setLevel(RDLogger.CRITICAL)

    def run(self):
        self.logger.info("Running populate monomers script")

        super().run()
        self.create_tables()
        self.insert_raw_smiles_data(use_cache=False)
        self.insert_ligand_counts(use_cache=False)
        self.insert_ligand_data()
        self.filter_ligands()
        self.fix_counts()

        self.logger.info("Finished populate monomers script")

    def create_tables(self):
        """
        Creates the tables that will store ligand definitions and counts
        """
        self.logger.info("Creating tables")

        self.run_commands("""
            create table {raw_data} (pdbid text, name text, formula text, smiles text, smilesoa text, inchi text);
            create table {counts}  (pdbid text, instances int default 0, deposits int default 0);
        """.format(raw_data=self.RAW_DATA_DB, counts=self.COUNTS_DB))

    def insert_raw_smiles_data(self, use_cache=False):
        """
        Downloads fresh smiles and inch data and inserts them into the database.
        """
        self.logger.info("Inserting raw smiles data")
        monomer_file = self._get_monomer_data(use_cache)
        self.copy_file_to_db(monomer_file, self.RAW_DATA_DB, columns="(pdbid, name, formula, smiles, smilesoa, inchi)")

    def _get_monomer_data(self, use_cache=False):
        smiles_source_url = "ftp://ftp.rcsb.org/pub/pdb/data/monomers/components.cif.gz"
        pdb_monomer_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "pdb_monomer_data.txt")
        result_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "monomers.tsv")

        if not use_cache or not os.path.exists(result_file):
            self.download_file(smiles_source_url, pdb_monomer_file, unpack=True)

            combined_df = self._parse_monomer_file(pdb_monomer_file)
            combined_df.to_csv(result_file, sep="\t")

        return result_file

    def _parse_monomer_file(self, pdb_monomer_file):
        with open(pdb_monomer_file, 'r+') as file:
            data = file.read()

            pdbids = [id.split('_')[1].strip() for id in re.compile("data_[A-Z0-9 ]{1,3}").findall(data)]
            names = [name.strip().strip('"').strip() for name in
                     re.compile("_chem_comp\.name[ ;\n]+(.+)", re.MULTILINE).findall(data)]
            formulas = [formula.strip().strip('"').strip() for formula in
                     re.compile("_chem_comp\.formula[ ;\n]+(.+)", re.MULTILINE).findall(data)]
            smiles_cactus = [[smile.split()[0], smile.split()[-1].strip('"').strip().lstrip(";")]
                             for smile in re.compile("[A-Z0-9 ]{1,3}[ ]+SMILES_CANONICAL[ ]+CACTVS[ ]+[0-9\.a-zA-Z]+[ ;\n]+.+",
                                                     re.MULTILINE).findall(data)]
            smiles_openeye = [[smile.split()[0], smile.split()[-1].strip('"').strip().lstrip(";")]
                             for smile in re.compile("[A-Z0-9 ]{1,3}[ ]+SMILES_CANONICAL[ ]+\"OpenEye OEToolkits\"[ ]+[0-9\.a-zA-Z]+[ ;\n]+.+",
                                                     re.MULTILINE).findall(data)]
            inchs = [[inch.split()[0], inch.split()[-1].strip('"').strip().lstrip(";")] for inch in
                     re.compile("[A-Z0-9 ]{1,3}[ ]+InChI[ ]+InChI[ ]+[0-9\.a-zA-Z]+[ ;\n]+.+", re.MULTILINE).findall(data)]
            release_statuses = [status.split()[1] for status in
                                re.compile("_chem_comp\.pdbx_release_status[ ]+[A-Z_0-9]+").findall(data)]

        names_df = pd.DataFrame(data={"pdbid": pdbids, "name": names, "formula": formulas})
        names_df = names_df.set_index('pdbid')

        cactus_df = pd.DataFrame(data=smiles_cactus, columns=["pdbid", "smiles"])
        cactus_df = cactus_df.set_index('pdbid')
        openeye_df = pd.DataFrame(data=smiles_openeye, columns=["pdbid", "smilesoa"])
        openeye_df = openeye_df.set_index('pdbid')
        smiles_df = cactus_df.join(openeye_df, how="left")

        inchs_df = pd.DataFrame(data=inchs, columns=["pdbid", "inchi"])
        inchs_df = inchs_df.set_index('pdbid')

        status_df = pd.DataFrame(data={"pdbid": pdbids, "status": release_statuses})
        status_df = status_df.set_index('pdbid')

        combined_df = names_df.join(smiles_df, how="left").join(inchs_df, how="left").join(status_df, how="left")
        combined_df = combined_df.loc[combined_df.status != 'OBS', :]
        combined_df = combined_df.drop("status", axis=1)

        return combined_df

    def insert_ligand_counts(self, use_cache=False):
        """
        Gets ligand deposit and instance counts and inserts them into the database.
        """
        self.logger.info("Inserting ligand counts")

        instance_counts_file = self._download_instance_counts_from_pdbj(use_cache)
        self.copy_file_to_db(instance_counts_file, self.COUNTS_DB, columns="(pdbid, deposits, instances)")

    def _download_instance_counts_from_pdbj(self, use_cache=False):
        instance_counts_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "instances.txt")

        if not use_cache or not os.path.exists(instance_counts_file):
            self.query_pdbj("""
                select np.comp_id, count(np.pdbid) deposits, cast(sum(e.pdbx_number_of_molecules) as integer) instances
                from pdbj.entity e, pdbj.pdbx_entity_nonpoly np
                where e.pdbid = np.pdbid and e.id = np.entity_id
                group by comp_id
                order by instances desc
            """, "tsv", instance_counts_file)

        return instance_counts_file

    def insert_ligand_data(self):
        """
        Import (parse/convert) smiles of ligands
        """
        self.logger.info("Inserting ligand data")

        self.run_commands("""
            select * into {ligands} from (select {raw_data}.pdbid, mol_from_smiles({raw_data}.smiles::cstring) molecule 
            from {raw_data}) tmp;
            
            update {ligands} set molecule = NULL where molecule::text like '%.%';
        
            alter table {ligands} add column atoms integer;
            alter table {ligands} add column rings integer;
            alter table {ligands} add column aromatic_rings integer;
            alter table {ligands} add column weight integer;
        
            create index ligands_molecule_idx on {ligands} using gist(molecule);
            create index ligands_properties_idx on {ligands}(atoms, rings, aromatic_rings);
        """.format(ligands=self.LIGANDS_DB, raw_data=self.RAW_DATA_DB, counts=self.COUNTS_DB))

        self.parse_problematic_ligands()

        self.run_commands("""
            update {ligands} set(atoms, rings, aromatic_rings, weight) = 
            (mol_numheavyatoms(molecule), mol_numrings(molecule), mol_numaromaticrings(molecule), mol_amw(molecule));
        """.format(ligands=self.LIGANDS_DB))

    def _update_parsed_molecules(self, parsed):
        for pdbid, molecule in parsed:
            self.cursor.execute("""
                 update {ligands} 
                 set molecule = mol_from_pkl(%s) 
                 where pdbid = %s;
            """.format(ligands=self.LIGANDS_DB), (bytearray(molecule), pdbid))
        self.conn.commit()

    def parse_problematic_ligands(self):
        self.logger.info("Attempting to re-parse problematic molecules")

        problematic_df = pd.read_sql_query("""
            select pdbid, smiles, smilesoa, inchi
            from {raw_data}
            where pdbid in (select pdbid from {ligands} where molecule is null);
        """.format(raw_data=self.RAW_DATA_DB, ligands=self.LIGANDS_DB), self.conn)

        skipped_smiles, fixed_smiles = self.parse_with_fixes(problematic_df, try_sanitized=False)
        skipped_smilesoa, fixed_smilesoa = self.parse_with_fixes(
            problematic_df[problematic_df.pdbid.isin(skipped_smiles)],
            parser=Chem.MolFromSmiles,
            field='smilesoa'
        )
        skipped_inchi, fixed_inchi = self.parse_with_fixes(
            problematic_df[problematic_df.pdbid.isin(skipped_smilesoa)],
            parser=Chem.MolFromInchi,
            field='inchi'
        )

        self._update_parsed_molecules(fixed_smiles+fixed_smilesoa+fixed_inchi)
        self.logger.info("%d ligands omitted due to problems with interpretation: %s", len(skipped_inchi), skipped_inchi)

    def parse_with_fixes(self, df, parser=Chem.MolFromSmiles, field='smiles', try_sanitized=True):
        skipped = []
        fixed = []
        h_problems = []
        self.logger.info("Using '%s' parser with '%s' field to apply fixes", parser.__name__, field)
        for row in df.itertuples():
            string_rep = str(getattr(row, field))
            molecule = None
            try:
                if try_sanitized:
                    molecule = parser(string_rep, sanitize=True)
                    if molecule is None:
                        self.logger.debug("%s could not be sanitized: %s", row.pdbid, string_rep)

                if molecule is None:
                    molecule = parser(string_rep, sanitize=False)

                if molecule is None:
                    self.logger.debug("%s could not be parsed: %s", row.pdbid, string_rep)
                    skipped.append(row.pdbid)
                    continue

                fragments = Chem.GetMolFrags(molecule)

                if len(fragments) > 2:
                    self.logger.debug("Skipped %s because it contains more than two fragments: %s",
                                      row.pdbid, string_rep)
                    skipped.append(row.pdbid)
                    continue

                elif len(fragments) == 2 and len(fragments[-1]) < 2:
                    self.logger.debug("Removing a single atom fragment %s: %s",
                                      row.pdbid, string_rep)
                    atomToRemove = list(fragments[-1])[0]
                    editable_molecule = Chem.EditableMol(molecule)
                    editable_molecule.RemoveAtom(atomToRemove)
                    molecule = editable_molecule.GetMol()
                    Chem.SanitizeMol(molecule)

                try:
                    molecule = Chem.RemoveHs(molecule)
                except:
                    self.logger.debug("Could not remove Hs from %s, trying to fix", row.pdbid)
                    try:
                        molecule = self.fix_ligand(row.pdbid, molecule)
                        if molecule is None:
                            continue
                    except Chem.MolSanitizeException as e:
                        self.logger.debug("Unable to fix %s: %s", row.pdbid, e)
                        h_problems.append(row.pdbid)
                        skipped.append(row.pdbid)
                        continue

                fixed.append((row.pdbid, molecule.ToBinary()))
                self.logger.debug("Successfully parsed %s", row.pdbid)
            except Exception as e:
                self.logger.error("Problem with %s: %s", row.pdbid, e)
                skipped.append(row.pdbid)
                continue

        self.logger.info("%d ligands omitted due to problems with removing H atoms: %s", len(h_problems), h_problems)

        return skipped, fixed

    def fix_ligand(self, pdbid, molecule):
        for fix_name, pattern, fix in LIGAND_FIXES:
            matches = molecule.GetSubstructMatches(pattern)
            for match in matches:
                molecule = fix(molecule, match)
                if molecule is None:
                    self.logger.debug("%s fix %s failed", pdbid, fix_name)
                    raise Chem.MolSanitizeException('Unable to apply fix %s' % fix_name)
                elif molecule is False:
                    self.logger.debug("%s skipped by %s", pdbid, fix_name)
                    return None
                else:
                    self.logger.debug("%s applied fix %s", pdbid, fix_name)

        molecule = Chem.RemoveHs(molecule)
        return molecule

    def filter_ligands(self):
        """
        Filters out amino_acid versions (AAA_XXX) and ambiguous residues (ASX, GLX)
        """
        self.logger.info("Removing ambiguous residues")

        self.run_commands("""
            delete from {ligands} where pdbid like '%\_%';
            delete from {ligands} where pdbid in ('ASX', 'GLX');
            -- These give a problems - they for some reason match even if matching the substructure is shorter.
            delete from {ligands} where pdbid in ('TSD', 'SPW', 'NWN', '8WV', 'CUV', '9RD', 'AUF', 'E9C', 'LN8');
        """.format(ligands=self.LIGANDS_DB))

    def fix_counts(self):
        self.logger.info("Setting counts for structures without free ligands in the PDB")

        self.run_commands("""
            insert into {counts}
                select l.pdbid, 0, 0 
                from ligands l left outer join counts c on c.pdbid = l.pdbid 
                where c.instances is null;
        """.format(counts=self.COUNTS_DB))


# Main script #
if __name__ == '__main__':
    script = PopulateLigandsScript()
    script.run()
