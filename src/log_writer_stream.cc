#include "log_writer_stream.hpp"

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/filtering_stream.hpp>

namespace {
  const std::string getDateString() {
    time_t timer = time(NULL);
    tm* date = localtime(&timer);
    std::string dateStr = (boost::format("[%02d/%02d %02d:%02d:%02d]") 
			   % (date->tm_mon+1) % date->tm_mday
			   % date->tm_hour % date->tm_min % date-> tm_sec).str();
    return dateStr;
  }
}

namespace logs{
  io::filtering_ostream lout;
  LogType info("INFO");
  LogType debug("DEBUG");
  LogType warn("WARN");
  LogType error("ERROR");
  void log_init(const std::string& filename, bool verbose) {
    static std::ofstream ofs(filename.c_str()); // log file
    static io::tee_filter<std::ostream>  fileFilt(ofs); // tee all passed data to logfile
    logs::lout.push(fileFilt); // 1st, tee off any data to the file (raw boost XML)
    if (verbose) {
      logs::lout.push(std::cout); // 2nd, tee off any data to cout (raw boost XML)
    }
  }
  void close() {
    io::close(logs::lout);
  }
  std::ostream& operator<<(std::ostream& os, const LogType& lt) {
    os << getDateString() << " [time=(" << time(NULL) << ")] [" << lt.type << "] ";
    return os;
  }
}
