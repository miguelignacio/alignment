#ifndef _LINEARLEASTSQUARESTESTCASE_H
#define _LINEARLEASTSQUARESTESTCASE_H

#include <TestCase.h>
#include <extensions/HelperMacros.h>

using namespace CppUnit;

class LinearLeastSquaresTestCase : public TestCase
{
  CPPUNIT_TEST_SUITE(LinearLeastSquaresTestCase);
  CPPUNIT_TEST(testLinearFit);
  CPPUNIT_TEST(testMatrixOperations);
  CPPUNIT_TEST(testChi2one);
  CPPUNIT_TEST(testChi2two);
  CPPUNIT_TEST_SUITE_END();

public:
  LinearLeastSquaresTestCase();
  virtual ~LinearLeastSquaresTestCase();
  void testLinearFit();
  void testMatrixOperations();
  void testChi2one();
  void testChi2two();
};

#endif
