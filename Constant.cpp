#include "Constant.h"

Constant::Constant() {
}

Constant::~Constant() {
}

int Constant::K = 10;
bool Constant::IS_WEIGHTED = false;
int Constant::COMMUNITY_POPULATION = 8;
const double Constant::PERCENTAGE_THRESHOLD = 0.5;
const double Constant::EPSILON = 0.2;
const double Constant::DELTA = 0.2;
const int Constant::NUM_THREAD = 4;
bool Constant::IS_BOUNDED_THRESHOLD = false;
const bool Constant::MODEL = false; // true : LT ; false : IC
const bool Constant::GCS = false;
const double Constant::EIG_TIME = 10800;
