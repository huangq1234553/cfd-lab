#ifndef SIM_UVCOMM_TEST_H
#define SIM_UVCOMM_TEST_H

const int U_FACTOR = 41;
const int V_FACTOR = 43;

int uvCommTest(int mpiRank, int mpiNumProc);
int uvCommTest23(int mpiRank, int mpiNumProc);
int uvCommTest32(int mpiRank, int mpiNumProc);

#endif //SIM_UVCOMM_TEST_H
