void runpropka_(int* numatoms, char* pdb, char* outname);

void runPKA(int numatoms, char* pdb, char* outname){
  runpropka_(&numatoms, pdb , outname);
}
