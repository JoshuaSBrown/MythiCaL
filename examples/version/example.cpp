#include <mythical/version.hpp>
#include <iostream>

using namespace mythical;
int main() {
  std::cout << "This file has been built with:" << std::endl;
  std::cout << PROJECT_NAME << " " << PROJECT_VER << std::endl;
  return 0;
}
