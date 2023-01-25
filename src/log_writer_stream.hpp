#include <iostream>
#include <fstream>
#include <time.h>
#include <boost/format.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#ifndef LOG_WRITER_STREAM_H_
#define LOG_WRITER_STREAM_H_

namespace io = boost::iostreams;

namespace logs{
  struct LogType {
    std::string type;
    LogType() {}
    LogType(const std::string& type) : type(type) {}
  };
  extern io::filtering_ostream lout;
  extern LogType info;
  extern LogType debug;
  extern LogType warn;
  extern LogType error;
  void log_init(const std::string& filename, bool verbose = true);
  void close();
  std::ostream& operator<<(std::ostream& os, const LogType& lt);
}

#endif
