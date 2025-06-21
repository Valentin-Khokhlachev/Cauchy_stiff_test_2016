#include "csvstream.h"



csvstream::csvstream(const std::string filename, const std::string separator)
    : fs_()
    , is_first_(true)
    , separator_(separator)
    , escape_seq_("\"")
    , special_chars_("\"")
{
    this->file_name = filename;
    fs_.exceptions(std::ios::failbit | std::ios::badbit);
    fs_.open(filename);
    qDebug() << '"' << this->file_name.c_str() << "\" opened!";
}

csvstream::~csvstream()
{
    flush();
    fs_.close();
    qDebug() << '"' << this->file_name.c_str() << "\" closed!";
}

void csvstream::flush()
{
    fs_.flush();
}

void csvstream::endrow()
{
    fs_ << std::endl;
    is_first_ = true;
}

csvstream &csvstream::operator <<(csvstream &(*val)(csvstream &))
{
    return val(*this);
}

csvstream& csvstream::operator << (const char * val)
{
    return write(escape(val));
}

csvstream& csvstream::operator << (const std::string & val)
{
    return write(escape(val));
}

template<typename T>
csvstream& csvstream::operator << (const T& val)
{
    return write(val);
}

template<typename T>
csvstream &csvstream::write(const T &val)
{
    if (!is_first_){
        fs_ << separator_;
    }else{
        is_first_ = false;
    }
    fs_ << val;
    return *this;
}

std::string csvstream::escape(const std::string &val)
{
    std::ostringstream result;
    // result << '"';
    std::string::size_type to, from = 0u, len = val.length();
    while ( (from < len) &&
            (std::string::npos != (to = val.find_first_of(special_chars_, from)) )){
        result << val.substr(from, to - from) << escape_seq_ << val[to];
        from = to + 1;
    }
    // result << val.substr(from) << '"';
    result << val.substr(from) ;
    return result.str();
}


