wget https://repo.anaconda.com/archive/Anaconda3-5.2.0-Linux-x86_64.sh
chmod +x Anaconda3-5.2.0-Linux-x86_64.sh
./Anaconda3-5.2.0-Linux-x86_64.sh
conda update -n base conda
conda create -n ligand_clustering
source activate ligand_clustering

conda install -c rdkit rdkit rdkit-postgresql
#conda install -c rdkit rdkit rdkit-postgresql95
initdb -D rdkit-postgresql
sed -i "s/#port = 5432/port = 55543/" rdkit-postgresql/postgresql.conf
pg_ctl -D rdkit-postgresql -l rdkit_postregs.log start
createdb pdbmonomers -p 55543

# left just iun case we have to revert back to rdkit-postgresql95
cd anaconda3/envs/rdkit/lib
#for i in ~/anaconda3/lib/libncursesw.*; do echo $i; done
ln -s ~/anaconda2/lib/libncursesw.so.6.0 .
ln -s ~/anaconda2/lib/libncursesw.so.6 .
ln -s ~/anaconda2/lib/libtinfow.so.6 .
ln -s ~/anaconda2/lib/libtinfow.so.6.0 .
ln -s ~/anaconda2/lib/libtinfow.so .
psql pdbmonomers -p 55543
#sudo apt install libxrender1
