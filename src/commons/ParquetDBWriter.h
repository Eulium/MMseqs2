#ifndef PARQUETWRITER_H
#define PARQUETWRITER_H
#include "MemoryTracker.h"
#include "Parameters.h"
#include "TranslateNucl.h"
#include "Orf.h"
#include <variant> //new
#include <carquet/include/carquet/carquet.h>

typedef std::pair<void*,carquet_physical_type> Record;
using CellValue = std::variant<std::string, double, int32_t, float>;

class ParquetDBWriter : public MemoryTracker{
public:
    ParquetDBWriter();
    void init(const char* resultFile,const std::vector<int>& outcodes);
    ~ParquetDBWriter();
    void checkStatus(const std::string message);
    void writeCell(const CellValue& record_value, const int record_num);
    void addSeqBasedOnAln(const int record_num, const char *seq, unsigned int offset,
                        const std::string &bt, bool reverse, bool isReverseStrand,
                        bool translateSequence, const TranslateNucl &translateNucl);
    void writeBatchToFile(const std::vector<int>& outcodes);
    void close();
private:
    carquet_status_t status;
    carquet_status_t writer_close_status = CARQUET_ERROR_NOT_OPEN;

    carquet_error_t err;
    carquet_schema_t* schema;
    carquet_writer_options_t opts;
    carquet_writer_t* writer;
    std::map<int, Record> record_batch;
    int record_counter;

    void* typToPointer(carquet_physical_type type);
    void outcodeToField(const int outcode);
    void writeColumn(int column_number);
    void freeRecordBatch();
    void clearRecordBatch();
};

#endif