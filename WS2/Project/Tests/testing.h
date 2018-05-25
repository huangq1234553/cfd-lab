#ifndef SIM_TESTING_H
#define SIM_TESTING_H

void assertEqual(double value, double expectation, char *testName);

int expectEqual(double value, double expectation, char *testName, char *errorMsgFmt, ...);

#endif //SIM_TESTING_H
