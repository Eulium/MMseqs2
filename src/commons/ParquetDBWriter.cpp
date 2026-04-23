#include "Debug.h"
#include "Util.h"
#include "Parameters.h"
#include "TranslateNucl.h"
#include "Orf.h"
#include "ParquetDBWriter.h"


//#include <variant> //new
#include <carquet/include/carquet/carquet.h>
#define ZSTD_STATIC_LINKING_ONLY
#include <zstd.h>


struct StringColumn {
    std::vector<std::string> storage;
    std::vector<carquet_byte_array_t> values;

    void add(const std::string& s) {
        storage.push_back(s);
    }

    void finalize(){
        values.reserve(storage.size());
        for (std::string& entry : storage){

            values.push_back(carquet_byte_array_t{reinterpret_cast<uint8_t*>(const_cast<char*>(entry.data())),static_cast<int32_t>(entry.size())});
        }
    }

    void clear() {
        values.clear();
        storage.clear();
    }
};


void ParquetDBWriter::checkStatus(const std::string message){
    if (status != carquet_status_t::CARQUET_OK){
        Debug(Debug::ERROR) << message << "\n";
        EXIT(EXIT_FAILURE);
    }
}

ParquetDBWriter::ParquetDBWriter(){

};

ParquetDBWriter::~ParquetDBWriter(){
    if (writer_close_status == CARQUET_ERROR_NOT_OPEN){
        // the instance was not initialised
        return;
    }
    freeRecordBatch();
    if (writer_close_status != CARQUET_ERROR_ALREADY_CLOSED){
        close();
    }
}

void ParquetDBWriter::init(const char* resultFile,const std::vector<int>& outcodes){
    err = CARQUET_ERROR_INIT;
    schema = carquet_schema_create(&err);
    record_counter = 0;
    
    if (!schema) {
        status = CARQUET_ERROR_INTERNAL;
        checkStatus("Schema initialisation error");
    }

    for (size_t i = 0; i < outcodes.size(); i++){
        // creates schema column and returns void pointer
        outcodeToField(outcodes[i]);
        checkStatus("Adding column to schema failed");
    }

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_LZ4_RAW;
    // opts.compression_level = 6;
    opts.row_group_size = 128 * 1024 * 1024;
    opts.page_size = 1024 * 1024;
    opts.write_statistics = true;
    opts.write_crc = true;
    opts.write_page_index = true;
    opts.write_bloom_filters = true;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN; // might want to change, depends on performance
    opts.created_by = "MMseqs parquet writer";

    writer = carquet_writer_create(resultFile, schema, &opts, &err);
    if (!writer) {
        carquet_writer_abort(writer);
        freeRecordBatch();
        status = CARQUET_ERROR_INTERNAL;
        checkStatus("Writer initialisation error");
    } else {
        writer_close_status = CARQUET_OK;
    }
    carquet_schema_free(schema); // writer copies schema, can be freed early

};

void* ParquetDBWriter::typToPointer(carquet_physical_type type){
    switch (type){
        case carquet_physical_type::CARQUET_PHYSICAL_BYTE_ARRAY:
            return reinterpret_cast<void*>(new StringColumn());
        case carquet_physical_type::CARQUET_PHYSICAL_DOUBLE:
            return reinterpret_cast<void*>(new std::vector<double>());
        case carquet_physical_type::CARQUET_PHYSICAL_INT32:
            return reinterpret_cast<void*>(new std::vector<int32_t>());
        case carquet_physical_type::CARQUET_PHYSICAL_FLOAT:
            return reinterpret_cast<void*>(new std::vector<float>());
        default:
            status = CARQUET_ERROR_INTERNAL;
            checkStatus("Create pointer for unkown field type");
            return nullptr;
    }
}

void ParquetDBWriter::outcodeToField(const int outcode){
    carquet_logical_type_t string_type = {};
    string_type.id = CARQUET_LOGICAL_STRING;
    switch (outcode) {
        case Parameters::OUTFMT_QUERY:{
            status = carquet_schema_add_column(schema, "query", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* query = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(query,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TARGET:{
            status =carquet_schema_add_column(schema, "target", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* target = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(target,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_EVALUE:{
            status =carquet_schema_add_column(schema, "evalue", CARQUET_PHYSICAL_DOUBLE, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* evalue = typToPointer(CARQUET_PHYSICAL_DOUBLE);
            record_batch[record_counter] = Record(evalue,CARQUET_PHYSICAL_DOUBLE);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_GAPOPEN:{
            status =carquet_schema_add_column(schema, "gapopen", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* gapopen = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(gapopen,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_FIDENT:{
            status =carquet_schema_add_column(schema, "fident", CARQUET_PHYSICAL_FLOAT, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* findent = typToPointer(CARQUET_PHYSICAL_FLOAT);
            record_batch[record_counter] = Record(findent,CARQUET_PHYSICAL_FLOAT);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_PIDENT:{
            status =carquet_schema_add_column(schema, "pident", CARQUET_PHYSICAL_FLOAT, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* pindent = typToPointer(CARQUET_PHYSICAL_FLOAT);
            record_batch[record_counter] = Record(pindent,CARQUET_PHYSICAL_FLOAT);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_NIDENT:{
            status =carquet_schema_add_column(schema, "nident", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* nident = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(nident,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QSTART:{
            status =carquet_schema_add_column(schema, "Qstart", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qstart = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(qstart,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QEND:{
            status =carquet_schema_add_column(schema, "qend", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qend = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(qend,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QLEN:{
            status =carquet_schema_add_column(schema, "qlen", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qlen = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(qlen,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TSTART:{
            status = carquet_schema_add_column(schema, "tstart", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* tstart = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(tstart,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TEND:{
            status = carquet_schema_add_column(schema, "tend", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* tend = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(tend,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TLEN:{
            status = carquet_schema_add_column(schema, "tlen", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* tlen = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(tlen,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_ALNLEN:{
            status = carquet_schema_add_column(schema, "alnlen", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* alnlen = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(alnlen,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_RAW:{
            status = carquet_schema_add_column(schema, "raw", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* raw = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(raw,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_BITS:{
            status = carquet_schema_add_column(schema, "score", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* score = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(score,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_CIGAR:{
            status = carquet_schema_add_column(schema, "cigar", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* cigar = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(cigar,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QSEQ:{
            status = carquet_schema_add_column(schema, "qseq", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qseq = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(qseq,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TSEQ:{
            status = carquet_schema_add_column(schema, "tseq", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* tseq = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(tseq,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QHEADER:{
            status = carquet_schema_add_column(schema, "qheader", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qheader = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(qheader,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_THEADER:{
            status = carquet_schema_add_column(schema, "theader", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* theader = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(theader,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QALN:{
            status = carquet_schema_add_column(schema, "qaln", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qaln = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(qaln,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TALN:{
            status = carquet_schema_add_column(schema, "taln", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* taln = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(taln,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_MISMATCH:{
            status = carquet_schema_add_column(schema, "mismatch", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* mismatch = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(mismatch,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QCOV:{
            status = carquet_schema_add_column(schema, "qcov", CARQUET_PHYSICAL_FLOAT, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qcov = typToPointer(CARQUET_PHYSICAL_FLOAT);
            record_batch[record_counter] = Record(qcov,CARQUET_PHYSICAL_FLOAT);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TCOV:{
            status = carquet_schema_add_column(schema, "tcov", CARQUET_PHYSICAL_FLOAT, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* tcov = typToPointer(CARQUET_PHYSICAL_FLOAT);
            record_batch[record_counter] = Record(tcov,CARQUET_PHYSICAL_FLOAT);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QSET:{
            status = carquet_schema_add_column(schema, "qset", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qset = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(qset,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QSETID:{
            status = carquet_schema_add_column(schema, "qsetid", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qsetid = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(qsetid,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TSET:{
            status = carquet_schema_add_column(schema, "tset", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qset = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(qset,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TSETID:{
            status = carquet_schema_add_column(schema, "tsetid", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* tsetid = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(tsetid,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TAXID:{
            status = carquet_schema_add_column(schema, "taxid", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* taxid = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(taxid,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TAXNAME:{
            status = carquet_schema_add_column(schema, "taxname", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* taxname = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(taxname,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TAXLIN:{
            status = carquet_schema_add_column(schema, "taxlin", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* taxlin = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(taxlin,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_EMPTY:{
            status = carquet_schema_add_column(schema, "empty", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* empty = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(empty,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QORFSTART:{
            status = carquet_schema_add_column(schema, "qorfstart", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qorfstart = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(qorfstart,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QORFEND:{
            status = carquet_schema_add_column(schema, "qorfend", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qorfend = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(qorfend,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TORFSTART:{
            status = carquet_schema_add_column(schema, "torfstart", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* torfstart = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(torfstart,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TORFEND:{
            status = carquet_schema_add_column(schema, "torfend", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* torfend = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(torfend,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_PPOS:{
            status = carquet_schema_add_column(schema, "ppos", CARQUET_PHYSICAL_FLOAT, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* ppos = typToPointer(CARQUET_PHYSICAL_FLOAT);
            record_batch[record_counter] = Record(ppos,CARQUET_PHYSICAL_FLOAT);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QFRAME:{
            status = carquet_schema_add_column(schema, "qframe", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qframe = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(qframe,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TFRAME:{
            status = carquet_schema_add_column(schema, "tframe", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* tframe = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(tframe,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        default:
            status = CARQUET_ERROR_INTERNAL;
            checkStatus("Tried to create schema for unkown field type, outcode number:" + std::to_string(outcode));
            break;
    }
}

void ParquetDBWriter::addSeqBasedOnAln(const int record_num, const char *seq, unsigned int offset,
                        const std::string &bt, bool reverse, bool isReverseStrand,
                        bool translateSequence, const TranslateNucl &translateNucl) {
    void* record_pointer = record_batch[record_num].first;
    StringColumn* carbyte_pointer = reinterpret_cast<StringColumn*>(record_pointer);
    std::string out; 
    unsigned int seqPos = 0;
    char codon[3];
    for (uint32_t i = 0; i < bt.size(); ++i) {
        char seqChar = (isReverseStrand == true) ? Orf::complement(seq[offset - seqPos]) : seq[offset + seqPos];
        if (translateSequence) {
            codon[0] = (isReverseStrand == true) ? Orf::complement(seq[offset - seqPos])     : seq[offset + seqPos];
            codon[1] = (isReverseStrand == true) ? Orf::complement(seq[offset - (seqPos+1)]) : seq[offset + (seqPos+1)];
            codon[2] = (isReverseStrand == true) ? Orf::complement(seq[offset - (seqPos+2)]) : seq[offset + (seqPos+2)];
            seqChar = translateNucl.translateSingleCodon(codon);
        }
        switch (bt[i]) {
            case 'M':
                out.append(1, seqChar);
                seqPos += (translateSequence) ?  3 : 1;
                break;
            case 'I':
                if (reverse == true) {
                    out.append(1, '-');
                } else {
                    out.append(1, seqChar);
                    seqPos += (translateSequence) ?  3 : 1;
                }
                break;
            case 'D':
                if (reverse == true) {
                    out.append(1, seqChar);
                    seqPos += (translateSequence) ?  3 : 1;
                } else {
                    out.append(1, '-');
                }
                break;
        }
    }
    carbyte_pointer->add(out);
}


void ParquetDBWriter::clearRecordBatch(){
    for (const auto& Record : record_batch) {
        void* record_pointer = Record.second.first;
        carquet_physical_type type = Record.second.second;
        switch (type){
            case CARQUET_PHYSICAL_BYTE_ARRAY:{
                StringColumn* carbyte_pointer_a = reinterpret_cast<StringColumn*>(record_pointer);
                carbyte_pointer_a->clear();
                break;
            }
            case CARQUET_PHYSICAL_DOUBLE:{
                std::vector<double>* carbyte_pointer_b = reinterpret_cast<std::vector<double>*>(record_pointer);
                carbyte_pointer_b->clear();
                break;
            }
            case CARQUET_PHYSICAL_INT32:{
                std::vector<int32_t>* carbyte_pointer_c = reinterpret_cast<std::vector<int32_t>*>(record_pointer);
                carbyte_pointer_c->clear();
                break;
            }
            case CARQUET_PHYSICAL_FLOAT:{
                std::vector<float>* carbyte_pointer_d = reinterpret_cast<std::vector<float>*>(record_pointer);
                carbyte_pointer_d->clear();
                break;
            }
            default:
                status = CARQUET_ERROR_INTERNAL;
                checkStatus("Clear record of unkown field type");
                break;
        }
    }
}

void ParquetDBWriter::freeRecordBatch(){
    for (const auto& record : record_batch) {
        void* record_pointer = record.second.first;
        carquet_physical_type type = record.second.second;
        switch (type){
            case CARQUET_PHYSICAL_BYTE_ARRAY:{
                StringColumn* carbyte_pointer_a = reinterpret_cast<StringColumn*>(record_pointer);
                delete carbyte_pointer_a; 
                break;
            }
            case CARQUET_PHYSICAL_DOUBLE:{
                std::vector<double>* carbyte_pointer_b = reinterpret_cast<std::vector<double>*>(record_pointer);
                delete carbyte_pointer_b;
                break;
            }
            case CARQUET_PHYSICAL_INT32:{
                std::vector<int32_t>* carbyte_pointer_c = reinterpret_cast<std::vector<int32_t>*>(record_pointer);
                delete carbyte_pointer_c;
                break;
            }
            case CARQUET_PHYSICAL_FLOAT:{
                std::vector<float>* carbyte_pointer_d = reinterpret_cast<std::vector<float>*>(record_pointer);
                delete carbyte_pointer_d;
                break;
            }
            default:
                status = CARQUET_ERROR_INTERNAL;
                checkStatus("Free record of unkown field type");
                break;
        }
    }
}


void ParquetDBWriter::writeCell(const std::string& record_value, const int record_num){
    void* record_pointer = record_batch[record_num].first;
    carquet_physical_type type = record_batch[record_num].second;
    StringColumn* carbyte_pointer_a = reinterpret_cast<StringColumn*>(record_pointer);
    const std::string& value = record_value;
    carbyte_pointer_a->add(value);
}

void ParquetDBWriter::writeCell(const double& record_value, const int record_num){
    void* record_pointer = record_batch[record_num].first;
    carquet_physical_type type = record_batch[record_num].second;
    std::vector<double>* carbyte_pointer_b = reinterpret_cast<std::vector<double>*>(record_pointer);
    carbyte_pointer_b->emplace_back(record_value);
}

void ParquetDBWriter::writeCell(const int& record_value, const int record_num){
    void* record_pointer = record_batch[record_num].first;
    carquet_physical_type type = record_batch[record_num].second;
    std::vector<int32_t>* carbyte_pointer_c = reinterpret_cast<std::vector<int32_t>*>(record_pointer);
    carbyte_pointer_c->emplace_back(record_value);
}

void ParquetDBWriter::writeCell(const float& record_value, const int record_num){
    void* record_pointer = record_batch[record_num].first;
    carquet_physical_type type = record_batch[record_num].second;
    std::vector<float>* carbyte_pointer_d = reinterpret_cast<std::vector<float>*>(record_pointer);
    carbyte_pointer_d->emplace_back(record_value);
}



void ParquetDBWriter::writeColumn(int column_number){
    
    void* record_pointer = record_batch[column_number].first;
    carquet_physical_type type = record_batch[column_number].second;
    
    switch (type){
        case CARQUET_PHYSICAL_BYTE_ARRAY: {
            StringColumn* carbyte_pointer_a = reinterpret_cast<StringColumn*>(record_pointer);
            carbyte_pointer_a->finalize(); 
            status = carquet_writer_write_batch(writer, column_number, carbyte_pointer_a->values.data(),static_cast<int64_t>(carbyte_pointer_a->values.size()), NULL, NULL);
            checkStatus("Column write" + std::to_string(column_number) + "error");
            break;
        }
        case CARQUET_PHYSICAL_DOUBLE:{
            std::vector<double>* carbyte_pointer_b = reinterpret_cast<std::vector<double>*>(record_pointer);
            status = carquet_writer_write_batch(writer, column_number, carbyte_pointer_b->data(),static_cast<int64_t>(carbyte_pointer_b->size()), NULL, NULL);
            checkStatus("Column write" + std::to_string(column_number) + "error");
            break;
        }
        case CARQUET_PHYSICAL_INT32:{
            std::vector<int32_t>* carbyte_pointer_c = reinterpret_cast<std::vector<int32_t>*>(record_pointer);
            status = carquet_writer_write_batch(writer, column_number, carbyte_pointer_c->data(),static_cast<int64_t>(carbyte_pointer_c->size()), NULL, NULL);
            checkStatus("Column write" + std::to_string(column_number) + "error");
            break;
        }
        case CARQUET_PHYSICAL_FLOAT:{
            std::vector<float>* carbyte_pointer_d = reinterpret_cast<std::vector<float>*>(record_pointer);
            status = carquet_writer_write_batch(writer, column_number, carbyte_pointer_d->data(),static_cast<int64_t>(carbyte_pointer_d->size()), NULL, NULL);
            checkStatus("Column write" + std::to_string(column_number) + "error");
            break;
        }
        default:
            status = carquet_status::CARQUET_ERROR_INTERNAL;;
            checkStatus("Column write" + std::to_string(column_number) + "error"); 
            break;
    }
}

void ParquetDBWriter::writeBatchToFile(const std::vector<int>& outcodes){
    for(size_t i = 0; i < outcodes.size(); i++){
        writeColumn(i);
    }
    clearRecordBatch();
}

void ParquetDBWriter::close(){
    if (writer_close_status == CARQUET_ERROR_ALREADY_CLOSED){
        Debug(Debug::WARNING) << "Writer already closed" << "\n";
        return;
    }
    writer_close_status = carquet_writer_close(writer);
    if (writer_close_status != carquet_status_t::CARQUET_OK){
        carquet_writer_abort(writer);
        writer_close_status = CARQUET_ERROR_ALREADY_CLOSED;
        checkStatus("Closing parquet writer failed");
        
    } else {
        writer_close_status = CARQUET_ERROR_ALREADY_CLOSED;
    }
}

