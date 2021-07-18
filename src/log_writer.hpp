
#include <iostream>

#ifndef LOG_WRITER_H_
#define LOG_WRITER_H_
namespace logger{
  struct logLevel{
    enum TypeEnum{
      ALL   = 0x00000000,
      DEBUG = 0x00000001,
      INFO  = 0x00000010,
      WARN  = 0x00000100,
      ERROR = 0x00001000
    };
  };

  void setLogFileName(std::string filename);
  void setLogLevel(int level);
  void writeLog(std::string msg, logger::logLevel::TypeEnum level);
}

#endif
