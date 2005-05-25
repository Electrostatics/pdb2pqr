void propka_(int* numatoms, int* strlen, char* pdb, char* outname);

char* runPKA(char* pdb, int strlen, int numatoms, char* outname){
  propka_(&numatoms, &strlen, pdb , outname);
  return pdb;
}
