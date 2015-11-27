#include <string>
#include <sstream>
#include <stdlib.h>     /* atof */

#include "ConfigDetails.h"

ConfigDetails::ConfigDetails(std::string const& value) {
  value_=value;
}

#include <iostream>

ConfigDetails::ConfigDetails(const char* c) {
  value_=c;
}

ConfigDetails::ConfigDetails(double d) {
  std::stringstream s;
  s<<d;
  value_=s.str();
}

ConfigDetails::ConfigDetails(ConfigDetails const& other) {
  value_=other.value_;
}

ConfigDetails& ConfigDetails::operator=(ConfigDetails const& other) {
  value_=other.value_;
  return *this;
}

ConfigDetails& ConfigDetails::operator=(double i) {
  std::stringstream s;
  s << i;
  value_ = s.str();
  return *this;
}

ConfigDetails& ConfigDetails::operator=(std::string const& s) {
  value_=s;
  return *this;
}

ConfigDetails::operator std::string() const {
  return value_;
}

ConfigDetails::operator double() const {
  return atof(value_.c_str());
}
