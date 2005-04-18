void propka_(int* numatoms, int* strlen, char* pdb);

char* runPKA(char* pdb, int strlen, int linelen){
  int numatoms;
  if ((strlen % linelen) != 0){
      return "Error Occurred";
  }
  numatoms = (int)(strlen/linelen);
  propka_(&numatoms, &strlen, pdb);
  return pdb;
}
