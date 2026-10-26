#ifndef POCKETS_INITIAL_STATE_H
#define POCKETS_INITIAL_STATE_H
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#if defined(__unix__) || defined(__APPLE__)
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

// Read-only source; the model still owns its separate evolved fields.
// PKIN's 20-byte header misaligns doubles, so memcpy avoids unaligned C++
// pointer dereferences. Mapping pages are reclaimable, not a heap copy.
class PocketInitialState {
    std::size_t count_=0;
#if defined(__unix__) || defined(__APPLE__)
    const char* mapping_=nullptr;
    std::size_t bytes_=0;
#else
    std::vector<double> values_;
#endif
public:
    PocketInitialState(const std::string& path,int n,double box):count_(4*std::size_t(n)*n*n){
        static_assert(std::endian::native==std::endian::little);
        char magic[8];std::uint32_t dimension=0;double length=0;
#if defined(__unix__) || defined(__APPLE__)
        const int fd=::open(path.c_str(),O_RDONLY);
        if(fd<0)throw std::runtime_error("Cannot open initial state");
        struct stat info{};const auto expected=20+count_*sizeof(double);
        if(::fstat(fd,&info)!=0||info.st_size<0||std::uintmax_t(info.st_size)!=expected){
            ::close(fd);throw std::invalid_argument("Initial state byte count");
        }
        void* mapped=::mmap(nullptr,expected,PROT_READ,MAP_PRIVATE,fd,0);::close(fd);
        if(mapped==MAP_FAILED)throw std::runtime_error("Cannot map initial state");
        mapping_=static_cast<const char*>(mapped);bytes_=expected;
        std::memcpy(magic,mapping_,8);std::memcpy(&dimension,mapping_+8,4);std::memcpy(&length,mapping_+12,8);
#else
        std::ifstream in(path,std::ios::binary);
        in.read(magic,8);in.read(reinterpret_cast<char*>(&dimension),4);in.read(reinterpret_cast<char*>(&length),8);
        if(!in)throw std::invalid_argument("Invalid initial state header");
        values_.resize(count_);in.read(reinterpret_cast<char*>(values_.data()),count_*8);
        if(!in||in.peek()!=std::char_traits<char>::eof())throw std::invalid_argument("Initial state byte count");
#endif
        try{
            if(std::string(magic,8)!="PKIN0001"||dimension!=std::uint32_t(n)||length!=box)
                throw std::invalid_argument("Invalid PKIN0001 initial state/grid");
            for(std::size_t i=0;i<count_;++i)if(!std::isfinite((*this)[i]))
                throw std::invalid_argument("Nonfinite initial state");
        }catch(...){release();throw;}
    }
    PocketInitialState(const PocketInitialState&)=delete;
    PocketInitialState& operator=(const PocketInitialState&)=delete;
    ~PocketInitialState(){release();}
    double operator[](std::size_t i) const {
#if defined(__unix__) || defined(__APPLE__)
        double value;std::memcpy(&value,mapping_+20+i*sizeof(double),sizeof(double));return value;
#else
        return values_[i];
#endif
    }
    void release() noexcept {
#if defined(__unix__) || defined(__APPLE__)
        if(mapping_){::munmap(const_cast<char*>(mapping_),bytes_);mapping_=nullptr;bytes_=0;}
#else
        std::vector<double>().swap(values_);
#endif
    }
};
#endif
