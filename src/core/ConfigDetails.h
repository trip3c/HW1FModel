#ifndef CONFIGDETAILS_H__
#define CONFIGDETAILS_H__

#include <string>

class ConfigDetails {
public:
  ConfigDetails() {};
  explicit ConfigDetails(const std::string&);
  explicit ConfigDetails(double);
  explicit ConfigDetails(const char*);

  ConfigDetails(const ConfigDetails&);
  ConfigDetails& operator=(ConfigDetails const&);

  ConfigDetails& operator=(double);
  ConfigDetails& operator=(std::string const&);

public:
  operator std::string() const;
  operator double     () const;
private:
  std::string value_;
};

#endif
