#include <iomanip>
#include <string>
#include <iostream>
#include <limits>
#include "json.hpp"

typedef std::numeric_limits< double > dbl;

using namespace std;

using json = nlohmann::json;

void secantTable(string);