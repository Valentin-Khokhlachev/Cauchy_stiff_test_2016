#ifndef CSVSTREAM_H
#define CSVSTREAM_H

#include <QDebug>

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

class csvstream;

inline static csvstream& endrow(csvstream& file);
inline static csvstream& flush(csvstream& file);

class csvstream
{
    std::ofstream fs_;
    bool is_first_;
    const std::string separator_;
    const std::string escape_seq_;
    const std::string special_chars_;

public:
    /*!
     * \brief csvstream - конструктор потока вывода в csv файл
     * \param filename - имя файла
     * \param separator - разделитель (по умолчанию это ";")
     */
    csvstream(const std::string filename, const std::string separator = ";");

    /*!
     * \brief ~csvstream - деструктор потока вывода в csv файл
     */
    ~csvstream();

    /*!
     * \brief flush - "промывка" потока вывода в csv файл
     */
    void flush();

    /*!
     * \brief endrow - конец строки в csv файле
     */
    void endrow();

    /*!
     * \brief operator << - перегрузка потока вывода (1)
     * \return
     */
    csvstream& operator << (csvstream& (* val)(csvstream&));

    /*!
     * \brief operator << - перегрузка потока вывода (2)
     * \param val - C-style строка
     * \return
     */
    csvstream& operator << (const char * val);

    /*!
     * \brief operator << - перегрузка потока вывода (3)
     * \param val - строка из стандартной библиотеки C++
     * \return
     */
    csvstream& operator << (const std::string & val);

    /*!
     * \brief operator << - перегрузка потока вывода (4)
     * \param val - представитель типа T,
     * для которого (пере-)определена операция "<<"
     * \return
     */
    template<typename T>
    csvstream& operator << (const T& val);

private:
    /*!
     * \brief file_name - имя файла
     */
    std::string file_name = "\0";

    /*!
     * \brief write - главная функция-принтер
     * \param val - представитель типа T,
     * для которого (пере-)определена операция "<<"
     * \return
     */
    template<typename T>
    csvstream& write (const T& val);

    /*!
     * \brief escape - функция-обработчик escape-последовательностей (\ n, \ t,...)
     * \param val - исходная строка
     * \return обработанная строка
     */
    std::string escape(const std::string & val);
};

/*!
 * \brief endrow - конец строки csv файла
 * \param file - имя файла
 * \return
 */
inline static csvstream& endrow(csvstream& file)
{
    file.endrow();
    return file;
}

/*!
 * \brief flush - "промывка" потока вывода в csv файл
 * \param file - имя файла
 * \return
 */
inline static csvstream& flush(csvstream& file)
{
    file.flush();
    return file;
}

#endif // CSVSTREAM_H
