#ifndef __CONFIG_FILE_H__
#define __CONFIG_FILE_H__

#include <string>
#include <map>

#include "ConfigDetails.h"

class ConfigFile {
  std::map<std::string,ConfigDetails> content_;

public:
  ConfigFile(std::string const& configFile);

  ConfigDetails const& Value(std::string const& section, std::string const& entry) const;

  ConfigDetails const& Value(std::string const& section, std::string const& entry, double value);
  ConfigDetails const& Value(std::string const& section, std::string const& entry, std::string const& value);
};

#endif
