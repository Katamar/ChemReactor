module load gcc/8.3.0
rm ChemReactor
echo "remove old binary"
g++ -g ChemReactor.cpp -o ChemReactor
./ChemReactor
