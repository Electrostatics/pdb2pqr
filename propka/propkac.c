void propka_(int* numatoms, int* strlen, char* pdb);

char* runPKA(char* pdb, int strlen, int numatoms){
  propka_(&numatoms, &strlen, pdb);
  return pdb;
}
