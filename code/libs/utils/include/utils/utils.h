#pragma once

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <utils/logger.h>

// possible loss of data in conversion between double and float
#pragma warning(disable : 4224)
// deprecated/unsafe functions such as fopen
#pragma warning(disable : 4996)

#define GET_STRING_FROM_ARGUMENT_LIST(fmt, pBuffer)              \
    {                                                            \
        va_list args;                                            \
        pBuffer = nullptr;                                       \
        int length = 1024;                                       \
        int result = -1;                                         \
        while (result == -1) {                                   \
            delete[] pBuffer;                                    \
            pBuffer = new char[length + 1];                      \
            memset(pBuffer, 0, length + 1);                      \
            va_start(args, fmt);                                 \
            result = std::vsnprintf(pBuffer, length, fmt, args); \
            va_end(args);                                        \
            if (result >= length) result = -1;                   \
            length *= 2;                                         \
        }                                                        \
    }

#define RELEASE_STRING_FROM_ARGUMENT_LIST(pBuffer) \
    {                                              \
        delete[] pBuffer;                          \
        pBuffer = nullptr;                         \
    }

/**
 * This method throws an error with a specified text and arguments
 */
inline void throwError(const char *fmt, ...) {
    char *pBuffer = nullptr;
    GET_STRING_FROM_ARGUMENT_LIST(fmt, pBuffer);
    std::cout << "Error Thrown: " << pBuffer << std::endl;
    throw pBuffer;
    RELEASE_STRING_FROM_ARGUMENT_LIST(pBuffer);

#ifdef WIN32
    __debugbreak();
#endif
}

typedef struct key_word {
    std::string keywordText;
    int keywordID;
} Keyword;

inline bool fileExists(const char *fName) {
    FILE *fp = fopen(fName, "r");
    if (fp) {
        fclose(fp);
        return true;
    }
    return false;
}

class BasicTextParser {
protected:
    //stateful information...

    //buffer that holds the content of the entire file - stored explicitly so that we can look ahead, etc
    char *buffer = nullptr;
    long int bufferLength = 0;
    //as we parse the contents of the file, offset keeps track of where we are
    long int offset = 0;
    //for easily interpretable error messages, keep track of the row and column of the character we're currently parsing
    int row = 1;
    int col = 1;

    bool verbose = true;
    std::string fName;

#define parserError(fmt, ...) if (verbose) Logger::print("parser error (%s: %i,%i): " fmt "\n", fName.c_str(), row, col, ##__VA_ARGS__)


    bool endOfTextStream(){
        return peek() == '\0';
    }

    // Peek at the n-th character away from the stream's current offset.
    char peek(int n = 0) {
        if (offset + n < bufferLength)
            return buffer[offset + n];
        else
            return '\0';
    }

    // Checks if the current character matches the given input
    bool is(char m) {
        return m == peek();
    }

    // Checks if the current character matches any one of the given characters
    bool isOneOf(const char *m) {
        const char* res = strchr(m, peek());
        return res != NULL && *res != '\0';
    }

    // Checks the next characters in the stream to see if they match the prefix in a caseless way
    bool startsWith(const char *prefix) {
        if (offset + strlen(prefix) >= bufferLength)
            return false;

        const char *start = buffer + offset;

        while (*prefix) {
            if (tolower(*prefix) != tolower(*start))
                return false;
            prefix++;
            start++;
        }
        return true;
    }

    // Advances the stream forward n characters
    void increment(int n = 1) {
        for (int i = 0; i < n; i++) {
            if (peek() == '\n') {
                row++;
                col = 0;
            } else {
                col++;
            }
            offset++;
        }
    }

    bool isWhiteSpace(char ch){
        return ch == ' ' || ch == '\t' || ch == '\r' || ch == '\v' || ch == '\n';
    }

    // Parse any whitespace on the current line of text
    void parseWhitespace() {
        while (isOneOf(" \r\t\v"))
            increment();
    }

    // parse all white space and new lines, the cursor will be set either
    // on end of stream or the first upcoming non-white-space symbol
    void goToNextNonWhiteSpaceCharacter() {
        while (isOneOf(" \r\t\v\n"))
            increment();
    }

    // parses the current line in its entirety - useful when the rest of
    // the current line can be safely ignored, e.g. it's a comment
    void parseEntireLine() {
        while (!endOfTextStream() && peek() != '\n')
            increment();
    }

    // Parse the given string (in a non-case sensitive way).
    bool parseKeyword(const char *string) {
        if (startsWith(string)) {
            increment(strlen(string));
        } else {
            parserError("expected '%s'", string);
            return false;
        }
        return true;
    }

    // Parse the given string (in a non-case sensitive way).
    // Many BVH files don't respect case sensitivity so parsing any keywords
    // in a non-case sensitive way seems safer.
    bool parseIfKeyword(const char *string) {
        if (startsWith(string))
            if (isWhiteSpace(peek(strlen(string)))){
                increment(strlen(string));
                return true;
            }
        return false;
    }

    // Parse any whitespace followed by a newline
    bool parseNewline() {
        parseWhitespace();

        if (is('\n')) {
            increment();
            parseWhitespace();
            return true;
        } else {
            parserError("expected newline");
            return false;
        }
    }

    // Parse any whitespace and then a  string/name made up of the characters in cList
    bool parseName(std::string &name, const char* cList = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_:-./" ) {
        parseWhitespace();

        name = "";

        while (isOneOf(cList)) {
            name += peek();
            increment();
        }

        if (name.length() > 0) {
            parseWhitespace();
            return true;
        } else {
            parserError("expected name");
            return false;
        }
    }

    // Parse a double
    bool parseDouble(double &out) {
        parseWhitespace();
        out = 0;

        char *end;
        errno = 0;
        out = strtod(buffer + offset, &end);

        if (errno == 0) {
            increment(end - (buffer + offset));
            return true;
        } else {
            parserError("expected number (double)");
            return false;
        }
    }

    // Parse an integer value
    bool parseInt(int &out) {
        parseWhitespace();

        out = 0;

        char *end;
        errno = 0;
        out = (int) strtol(buffer + offset, &end, 10);

        if (errno == 0) {
            increment(end - (buffer + offset));
            return true;
        } else {
            parserError("expected number (int)");
            return false;
        }
    }

    //returns true if the load is succesful, false otherwise...
    bool loadBufferFromFile(const char *filename){
        this->fName = filename;
        freeBuffer();

        // Read file Contents
        FILE *f = fopen(filename, "rb");

        if (f == NULL) {
            if (verbose)
                Logger::print("Error: Could not find file '%s'\n", filename);
            return false;
        }

        fseek(f, 0, SEEK_END);
        bufferLength = ftell(f);
        fseek(f, 0, SEEK_SET);
        buffer = (char *) malloc(bufferLength + 1);
        fread(buffer, 1, bufferLength, f);
        buffer[bufferLength] = '\n';
        fclose(f);

        return true;
    }

    void freeBuffer(){
        free(buffer);
        buffer = nullptr;
    }

};

inline void exitWithError(const std::string& s){
    std::cout << s;
    exit(0);
}

