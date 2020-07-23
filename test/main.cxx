#include <stdlib.h>

#include <extensions/TestFactoryRegistry.h>
#include <extensions/HelperMacros.h>
#include <ui/text/TestRunner.h>

using namespace std;
using namespace CppUnit;

int main(int argc, char * argv[])
{
  TextTestRunner runner;
  TestFactoryRegistry & registry = TestFactoryRegistry::getRegistry();
  runner.addTest(registry.makeTest());
  runner.run();

  return EXIT_SUCCESS;
}
