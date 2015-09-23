#ifndef __PDF_EXCEPTIONS__
#define __PDF_EXCEPTIONS__

class PdfException : public std::runtime_error {
 public:
 PdfException(const std::string& errorStr) : runtime_error(errorStr) {}
};


class DimensionError : public PdfException {
 public:
 DimensionError(const std::string& errorStr) : PdfException(errorStr){}
};



class OutOfBoundsError : public PdfException{
 public:
 OutOfBoundsError(const std::string& errorStr) : PdfException(errorStr){}
};

class BinError : public PdfException{
 public:
 BinError(const std::string& errorStr): PdfException(errorStr) {}
};

#endif