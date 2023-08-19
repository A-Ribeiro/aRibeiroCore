#include "StringUtil.h"

#include <stdio.h>  // printf
#include <stdarg.h> //va_start
#include <stdlib.h>
#include <string.h>
#include <wchar.h> // wprintf

#include <algorithm> // transform

#include <cctype> // tolower

namespace aRibeiro
{

    // private copy constructores, to avoid copy...
    StringUtil::StringUtil(const StringUtil &v) {}
    void StringUtil::operator=(const StringUtil &v) {}

    StringUtil::StringUtil() {}

    const char *StringUtil::char_ptr() const
    {
        if (char_buffer.size() == 0)
            return NULL;
        return &char_buffer[0];
    }

    const wchar_t *StringUtil::wchar_ptr() const
    {
        if (wchar_buffer.size() == 0)
            return NULL;
        return &wchar_buffer[0];
    }

    const char *StringUtil::printf(const char *format, ...)
    {

        va_list args;

        va_start(args, format);
        char_buffer.resize(vsnprintf(NULL, 0, format, args) + 1);
        va_end(args);

        va_start(args, format);
        int len = vsnprintf(&char_buffer[0], char_buffer.size(), format, args);
        va_end(args);

        return char_ptr();
    }

    const wchar_t *StringUtil::wprintf(const wchar_t *format, ...)
    {
        va_list args;

        va_start(args, format);
        wchar_buffer.resize(vswprintf(NULL, 0, format, args) + 1);
        va_end(args);

        va_start(args, format);
        int len = vswprintf(&wchar_buffer[0], wchar_buffer.size(), format, args);
        va_end(args);

        return wchar_ptr();
    }

    std::wstring StringUtil::toWString(const std::string &str)
    {
        return std::wstring(str.begin(), str.end());
    }

    std::string StringUtil::toString(const std::wstring &wstr)
    {
        return std::string(wstr.begin(), wstr.end());
    }

    bool StringUtil::startsWith(const std::string str, const std::string prefix)
    {
        return ((prefix.size() <= str.size()) &&
                std::equal(prefix.begin(), prefix.end(), str.begin()));
    }

    bool StringUtil::endsWith(const std::string str, const std::string sufix)
    {
        return ((sufix.size() <= str.size()) &&
                std::equal(sufix.begin(), sufix.end(), str.begin() + (str.size() - sufix.size())));
    }

    std::string StringUtil::toLower(const std::string str)
    {
        std::string aux = str;
        std::transform(aux.begin(), aux.end(), aux.begin(), ::tolower);
        return aux;
    }

    std::string StringUtil::toUpper(const std::string str)
    {
        std::string aux = str;
        std::transform(aux.begin(), aux.end(), aux.begin(), ::toupper);
        return aux;
    }

    bool StringUtil::contains(const std::string str, const std::string v)
    {
        return strstr(str.c_str(), v.c_str()) != 0;
    }

    std::vector<std::string> StringUtil::tokenizer(const std::string &input, const std::string &delimiter)
    {

        int start = 0;

        const uint8_t *input_str = (const uint8_t *)input.c_str();
        int input_size = (int)input.length();
        const uint8_t *delimiter_str = (const uint8_t *)delimiter.c_str();
        int delimiter_size = (int)delimiter.length();

        std::vector<std::string> result;
        std::string result_str;
        while (start < input_size)
        {
            int index_of_delimiter = array_index_of(input_str, start, input_size, delimiter_str, delimiter_size);
            if (index_of_delimiter == start)
            {
                start += delimiter_size;
                result.push_back("");
                continue;
            }
            int delta = index_of_delimiter - start;
            result_str.resize(delta);
            // for (int i = start, count = 0; i < index_of_delimiter; i++, count++)
            //     result_str[count] = (char)input_str[i];
            memcpy(&result_str[0], &input_str[start], delta * sizeof(uint8_t));
            result.push_back(result_str);
            start = index_of_delimiter + delimiter_size;
        }

        if (result.size() == 0)
            result.push_back("");

        return result;
    }

    void StringUtil::replaceAll(std::string *data, const std::string &toSearch, const std::string &replaceStr)
    {
        // Get the first occurrence
        size_t pos = data->find(toSearch);
        // Repeat till end is reached
        while (pos != std::string::npos)
        {
            // Replace this occurrence of Sub String
            data->replace(pos, toSearch.size(), replaceStr);
            // Get the next occurrence from the current position
            pos = data->find(toSearch, pos + replaceStr.size());
        }
    }
    std::string StringUtil::trim(const std::string &str)
    {
        // return std::regex_replace(str, std::regex("(^[ \n]+)|([ \n]+$)"),"");
        std::string::const_iterator it = str.begin();
        while (it != str.end() && std::isspace(*it))
            it++;

        std::string::const_reverse_iterator rit = str.rbegin();
        while (rit.base() != it && std::isspace(*rit))
            rit++;

        return std::string(it, rit.base());
    }
    std::string StringUtil::quote_cmd(const std::string &str)
    {
        // return std::string("\"") + std::regex_replace(str, std::regex("[\\\\\"]"), "\\$&") + "\"";
        std::string result = str;
        StringUtil::replaceAll(&result, "\\", "\\\\");
        StringUtil::replaceAll(&result, "\"", "\\\"");
        result = std::string("\"") + result + "\"";
        return result;
    }

    std::vector<std::string> StringUtil::parseArgv(const std::string &_cmd)
    {
        std::string cmd = trim(_cmd);

        std::vector<std::string> result;
        enum states
        {
            normal,
            enter_string,
            enter_string_concatenate,
            include_next_str,
            normal_slash
        };
        states state = normal;
        std::string aux = "";
        for (int i = 0; i < cmd.size(); i++)
        {
            char chr = cmd[i];
            char next_chr = cmd[i];
            if (i < (int)cmd.size() - 1)
                next_chr = cmd[i + 1];

            if (state == normal_slash)
            {
                state = normal;
                aux += chr;
            }
            else if (state == normal && chr == '\\')
            {
                state = normal_slash;
                aux += chr;
            }
            else if (state == normal && chr == '"')
            {
                // string opening...
                aux = trim(aux);
                if (aux.size() > 0){
                    aux += chr;
                    state = enter_string_concatenate;
                }else{
                    //result.push_back(aux);
                    state = enter_string;
                    aux = "";
                }
            }
            else if (state == normal && chr == ' ')
            {
                // one new command
                aux = trim(aux);
                if (aux.size() > 0)
                    result.push_back(aux);
                aux = "";
            }
            else if (state == enter_string && (chr == '\\' && next_chr == '"' || chr == '\\' && next_chr == '\\'))
            {
                // skip slash for inner string
                state = include_next_str;
            }
            else if (state == include_next_str)
            {
                // back to enter string
                state = enter_string;
                aux += chr;
            }
            else if (state == enter_string && chr == '"')
            {
                // exit string
                state = normal;
                result.push_back(aux);
                aux = "";
            }
            else if (state == enter_string_concatenate && chr == '"')
            {
                // exit string concatenate
                aux += chr;
                state = normal;
                aux = trim(aux);
                result.push_back(aux);
                aux = "";
            }
            else if (state == enter_string || state == enter_string_concatenate || state == normal)
            {
                aux += chr;
            }
        }

        aux = trim(aux);
        if (aux.size() > 0)
        {
            result.push_back(aux);
            aux = "";
        }

        // ::printf("[StringUtil]Parsing ARGV:\n");
        // for(size_t i=0;i<result.size();i++){
        //     ::printf("    %i -> '%s'\n", (int)i, result[i].c_str());
        // }

        return result;
    }

    std::string StringUtil::argvToCmd(const std::vector<std::string> &argv)
    {
        std::string result;

        for (size_t i = 0; i < argv.size(); i++)
        {
            std::string cmd = trim(argv[i]);
            if (cmd.size() == 0)
            {
                if (result.size() > 0)
                    result += " ";
                result += "\"\"";
                // continue;
            }
            else if (cmd.find(" ") != std::string::npos)
            {
                // found space, need to quote the string
                if (result.size() > 0)
                    result += " ";
                result += quote_cmd(cmd);
            }
            else
            {
                if (result.size() > 0)
                    result += " ";
                result += cmd;
            }
        }

        return result;
    }

}
