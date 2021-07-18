#include <iostream>
#include <fstream>
#include <time.h>
#include <boost/format.hpp>
#include "log_writer.hpp"


namespace {
  std::string log_file_name = "";
  int log_write_level = 0;
  std::ofstream log_writer;

  void write(std::string msg, logger::logLevel::TypeEnum level, std::string filename){
    std::string levelStr;
    switch(level){
    case logger::logLevel::DEBUG:
      levelStr = "[DEBUG]";
      break;
    case logger::logLevel::INFO:
      levelStr = "[ INFO]";
      break;
    case logger::logLevel::WARN:
      levelStr = "[ WARN]";
      break;
    case logger::logLevel::ERROR:
      levelStr = "[ERROR]";
      break;
    case logger::logLevel::ALL:
      break;
    }

#pragma omp critical(logWrite)
    {
      log_writer.open(filename.c_str(), std::ios::app);
      if(!log_writer){
	std::cerr << "Can't open file: " << filename << std::endl;
      }else{
	
	time_t timer;
	tm* date;	
	timer = time(NULL);
	date = localtime(&timer);
	std::string dateStr = (boost::format("[%02d/%02d %02d:%02d:%02d]") 
			       % date->tm_mon % date->tm_mday
			       % date->tm_hour % date->tm_min % date-> tm_sec).str(); 
	log_writer << dateStr << levelStr << " " << msg << std::endl;

      }
      log_writer.close();
    }//end of critical section.
  }
}//namespace

namespace logger{

  void setLogFileName(std::string filename){
    // change once only.
    if(log_file_name.compare("") != 0){
      std::cerr << "Log file name is already set." << std::endl;
    }else{
      log_file_name = filename;
      log_writer.open(filename.c_str());
      log_writer.close();
      /*      if(!log_writer){
	std::cerr << "Can't open file: " << filename << std::endl;
      }else{
      }*/
      
    }
  }
  
  void setLogLevel(int level){
    if(log_write_level != 0){
      std::cerr << "Log level is already set." << std::endl;
    }else{
      log_write_level = level;
    }
  }

  void writeLog(std::string msg, logger::logLevel::TypeEnum level){
    if(log_write_level <= level){ 
      write(msg, level, log_file_name);
    }
  }
  
}//namespace log

/*
int main(void){
  log::setLogFileName("test.log");
  log::setLogLevel(log::logLevel::ALL);

  log::writeLog("test", log::logLevel::WARN);
}
*/
