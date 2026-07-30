/*
 * Copyright (C) 2011-2014 by
 * Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
 * ... and other members of the Reconstruction Team of David Haussler's
 * lab (BME Dept. UCSC).
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
#include <assert.h>
#include <errno.h> // file existence via ENOENT
#include <getopt.h>
#include <inttypes.h>
#include <limits.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include "common.h"
#include "sharedMaf.h"
#include "buildVersion.h"

const char *g_version = "version 0.1 September 2012";

void version(void);
void usage(void);
void parseOptions(int argc, char **argv, char *filename, char **nameList,
                  bool *isInclude, int64_t *blockDegLT, int64_t *blockDegGT,
                  double *maxRefNFrac);
void checkRegion(unsigned lineno, char *fullname, uint64_t pos, uint64_t start,
                 uint64_t length, uint64_t sourceLength, char strand);
bool nameOnList(char *name, char **namelist, unsigned n);
void reportBlock(mafBlock_t *mb, char **names, unsigned n, bool isInclude);
void checkBlock(mafBlock_t *mb, char **names, unsigned n, bool isInclude,
                int64_t excludeBlockDegreeGT, int64_t excludeBlockDegreeLT,
                double maxRefNFrac);
void filterInput(mafFileApi_t *mfa, char **names, unsigned n,
                 bool isInclude, int64_t excludeBlockDegreeGT,
                 int64_t excludeBlockDegreeLT, double maxRefNFrac);
unsigned countNames(char *s);
char* readNameListFile(const char *filename);
char** extractNames(char *nameList, unsigned n);
void destroyNameList(char **names, unsigned n);
void reportNames(char **names, unsigned n);

void version(void) {
    fprintf(stderr, "mafFilter, %s\nbuild: %s, %s, %s\n\n", g_version, g_build_date,
            g_build_git_branch, g_build_git_sha);
}
void usage(void) {
    version();
    fprintf(stderr, "Usage: mafFilter --maf [path to maf] "
            "[options]\n\n"
            "mafFilter is a program that will look through a\n"
            "maf file block by block and excise out sequence lines\n"
            "that match criteria established by the user on the\n"
            "command line. For example one can filter out all\n"
            "sequence lines that start with 'hg18' using --excludeSeq\n"
            "or filter for sequence lines starting with only 'hg19',\n"
            "'mm9' and 'rn4' using --includeSeq.\n\n"
            "For long lists, --includeSeqFile and --excludeSeqFile read the\n"
            "same names from a file, one per line, which avoids having to\n"
            "build a comma separated list on the command line. Blank lines\n"
            "and lines beginning with '#' are ignored.\n\n");
    fprintf(stderr, "Options: \n");
    usageMessage('h', "help", "show this help message and exit.");
    usageMessage('m', "maf", "path to maf file. use - for stdin.");
    usageMessage('i', "includeSeq", "comma separated list of sequence names to include.");
    usageMessage('e', "excludeSeq", "comma separated list of sequence names to exclude.");
    usageMessage('I', "includeSeqFile", "file of sequence names to include, one per line.");
    usageMessage('E', "excludeSeqFile", "file of sequence names to exclude, one per line.");
    usageMessage('g', "noDegreeGT", "filter out all blocks with degree greater than this value.");
    usageMessage('l', "noDegreeLT", "filter out all blocks with degree less than this value.");
    usageMessage('N', "maxRefNFrac", "filter out all blocks whose reference sequence is more than this fraction of Ns");    
    usageMessage('v', "verbose", "turns on verbose output.");
    exit(EXIT_FAILURE);
}
void parseOptions(int argc, char **argv, char *filename, char **nameList, bool *isInclude, int64_t *blockDegGt, int64_t *blockDegLt, double *maxRefNFrac) {
    extern int g_debug_flag;
    extern int g_verbose_flag;
    int c;
    bool setMafName = false, setNames = false, setBlockLimits = false, setRefGap = false;
    while (1) {
        static struct option longOptions[] = {
            {"debug", no_argument, &g_debug_flag, 1},
            {"verbose", no_argument, 0, 'v'},
            {"help", no_argument, 0, 'h'},
            {"version", no_argument, 0, 0},
            {"maf",  required_argument, 0, 'm'},
            {"includeSeq",  required_argument, 0, 'i'},
            {"excludeSeq",  required_argument, 0, 'e'},
            {"includeSeqFile",  required_argument, 0, 'I'},
            {"excludeSeqFile",  required_argument, 0, 'E'},
            {"noDegreeGT", required_argument, 0, 'g'},
            {"noDegreeLT", required_argument, 0, 'l'},
            {"maxRefNFrac", required_argument, 0, 'N'},
            {0, 0, 0, 0}
        };
        int longIndex = 0;
        c = getopt_long(argc, argv, "m:i:e:I:E:g:l:v:N:h",
                        longOptions, &longIndex);
        if (c == -1) {
            break;
        }
        switch (c) {
        case 0:
            if (strcmp("version", longOptions[longIndex].name) == 0) {
                version();
                exit(EXIT_SUCCESS);
            }
            break;
        case 'm':
            setMafName = true;
            // bounded copy: an over-long path used to overrun filename
            if (snprintf(filename, kMaxStringLength, "%s", optarg) >= kMaxStringLength) {
                fprintf(stderr, "--maf path is longer than the %d character maximum\n", kMaxStringLength - 1);
                exit(EXIT_FAILURE);
            }
            break;
        case 'i':
            setNames = true;
            *isInclude = true;
            // the name list is duplicated onto the heap rather than read into a fixed buffer, as
            // a list naming every sequence of a large alignment runs to many kilobytes
            free(*nameList);
            *nameList = de_strdup(optarg);
            break;
        case 'e':
            setNames = true;
            *isInclude = false;
            free(*nameList);
            *nameList = de_strdup(optarg);
            break;
        case 'I':
            setNames = true;
            *isInclude = true;
            free(*nameList);
            *nameList = readNameListFile(optarg);
            break;
        case 'E':
            setNames = true;
            *isInclude = false;
            free(*nameList);
            *nameList = readNameListFile(optarg);
            break;
        case 'g':
            setNames = true;
            sscanf(optarg, "%" PRIi64, blockDegGt);
            break;
        case 'l':
            setNames = true;
            sscanf(optarg, "%" PRIi64, blockDegLt);
            break;
        case 'v':
            g_verbose_flag++;
            break;
        case 'N':
            setRefGap = true;
            sscanf(optarg, "%lf", maxRefNFrac);
            break;
        case 'h':
        case '?':
            usage();
            break;
        default:
            abort();
        }
    }
    if (!setMafName) {
        fprintf(stderr, "specify --maf\n");
        usage();
    }
    if ((!setRefGap && ((setNames && setBlockLimits) || !(setNames || setBlockLimits))) ||
        (setRefGap && (setNames || setBlockLimits))) {
        fprintf(stderr, "specify *one* from --maxRefNFrac or [--includeSeq or --excludeSeq] or [--noDegreeGT or --noDegreeLT]\n");
        usage();
    }
    // Check there's nothing left over on the command line
    if (optind < argc) {
        char errorString[30] = "Unexpected arguments:";
        while (optind < argc) {
            strcat(errorString, " ");
            strcat(errorString, argv[optind++]);
        }
        fprintf(stderr, "%s\n", errorString);
        usage();
    }
}
bool nameOnList(char *name, char **namelist, unsigned n) {
    for (unsigned i = 0; i < n; ++i) {
        if ((strncmp(name, namelist[i], strlen(namelist[i])) == 0)) {
            return true;
        }
    }
    return false;
}
void reportBlock(mafBlock_t *mb, char **names, unsigned n, bool isInclude) {
    // report the block being mindful of only including or excluding.
    mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
    while (ml != NULL) {
        if (maf_mafLine_getType(ml) != 's') {
            // report all sequence lines
            printf("%s\n", maf_mafLine_getLine(ml));
            ml = maf_mafLine_getNext(ml);
            continue;
        }
        if (n > 0) {
            if (isInclude) {
                if (nameOnList(maf_mafLine_getSpecies(ml), names, n)) {
                    printf("%s\n", maf_mafLine_getLine(ml));
                    ml = maf_mafLine_getNext(ml);
                    continue;
                }
            } else {
                if (!nameOnList(maf_mafLine_getSpecies(ml), names, n)) {
                    printf("%s\n", maf_mafLine_getLine(ml));
                    ml = maf_mafLine_getNext(ml);
                    continue;
                }
            }
        } else {
            // report entire block, this came from one of the blockDegree options
            printf("%s\n", maf_mafLine_getLine(ml));
        }
        ml = maf_mafLine_getNext(ml);
    }
    printf("\n");
}
void checkBlock(mafBlock_t *mb, char **names, unsigned n, bool isInclude,
                int64_t excludeBlockDegreeGT, int64_t excludeBlockDegreeLT,
                double maxRefNFrac) {
    // walk through the maf lines and see if this block should be reported
    mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
    while (ml != NULL) {
        if (maf_mafLine_getType(ml) != 's') {
            ml = maf_mafLine_getNext(ml);
            continue;
        }
        if (maxRefNFrac >= 0.) {
            char *sequence = maf_mafLine_getSequence(ml);
            int64_t length = strlen(sequence);
            int64_t bases = 0;
            int64_t Ns = 0;
            for (int64_t i = 0; i < length; ++i) {
                if (sequence[i] == 'N' || sequence[i] == 'n') {
                    ++Ns;
                }
                if (sequence[i] != '-') {
                    ++bases;
                }
            }
            if ((double)Ns / (double)bases <= maxRefNFrac) {
                reportBlock(mb, names, n, isInclude);
            } 
            return;
        } else if (n > 0) {
            // filtering on names
            if (isInclude) {
                if (nameOnList(maf_mafLine_getSpecies(ml), names, n)) {
                    reportBlock(mb, names, n, isInclude);
                    return;
                }
            } else {
                if (!nameOnList(maf_mafLine_getSpecies(ml), names, n)) {
                    reportBlock(mb, names, n, isInclude);
                    return;
                }
            }
        } else {
            // filtering on block degrees
            int64_t m = maf_mafBlock_getNumberOfSequences(mb);
            if (excludeBlockDegreeGT != -1 && excludeBlockDegreeLT != -1) {
                if (m >= excludeBlockDegreeLT && m <= excludeBlockDegreeGT) {
                    reportBlock(mb, names, n, isInclude);
                    return;
                }
            } else if (excludeBlockDegreeGT != -1) {
                if (m <= excludeBlockDegreeGT) {
                    reportBlock(mb, names, n, isInclude);
                    return;
                }
            } else {
                if (m >= excludeBlockDegreeLT) {
                    reportBlock(mb, names, n, isInclude);
                    return;
                }
            }
        }
        ml = maf_mafLine_getNext(ml);
    }
}
void filterInput(mafFileApi_t *mfa, char **names, unsigned n,
                 bool isInclude, int64_t excludeBlockDegreeGT,
                 int64_t excludeBlockDegreeLT, double maxRefNFrac) {
    mafBlock_t *thisBlock = NULL;
    bool headBlock = true;
    while ((thisBlock = maf_readBlock(mfa)) != NULL) {
        if (headBlock) {
            reportBlock(thisBlock, names, n, isInclude);
            headBlock = false;
            maf_destroyMafBlockList(thisBlock);
            continue;
        }
        checkBlock(thisBlock, names, n, isInclude, excludeBlockDegreeGT, excludeBlockDegreeLT, maxRefNFrac);
        maf_destroyMafBlockList(thisBlock);
    }
}
char* readNameListFile(const char *filename) {
    // Read sequence names one per line and return them as the comma separated string that
    // --includeSeq / --excludeSeq would have been given, so the rest of the program is unchanged.
    // Blank lines and lines starting with '#' are skipped, which lets a list carry a header or
    // be commented, and leading/trailing whitespace is trimmed.
    FILE *fh = de_fopen(filename, "r");
    int64_t lineLen = kMaxStringLength;
    char *line = (char*) de_malloc(lineLen);
    size_t capacity = kMaxStringLength, used = 0;
    char *list = (char*) de_malloc(capacity);
    list[0] = '\0';
    int64_t status = 0;
    while (status != -1) {
        status = de_getline(&line, &lineLen, fh);
        // de_getline signals end of file with -1 but still returns a final line that had no
        // trailing newline, so the buffer is always processed before the loop exits
        char *name = line;
        while (*name == ' ' || *name == '\t') {
            ++name;
        }
        char *end = name + strlen(name);
        while (end > name && (end[-1] == ' ' || end[-1] == '\t')) {
            --end;
        }
        *end = '\0';
        if (*name == '\0' || *name == '#') {
            continue;
        }
        size_t needed = used + strlen(name) + 2; // the name, a separating comma and the terminator
        if (needed > capacity) {
            while (needed > capacity) {
                capacity *= 2;
            }
            char *bigger = (char*) realloc(list, capacity);
            if (bigger == NULL) {
                fprintf(stderr, "(de_) realloc failed on a request for %zu bytes\n", capacity);
                exit(EXIT_FAILURE);
            }
            list = bigger;
        }
        if (used > 0) {
            list[used++] = ',';
        }
        strcpy(list + used, name);
        used += strlen(name);
    }
    free(line);
    fclose(fh);
    if (used == 0) {
        fprintf(stderr, "ERROR, %s does not contain any sequence names\n", filename);
        exit(EXIT_FAILURE);
    }
    return list;
}
unsigned countNames(char *s) {
    unsigned i, n;
    if (s == NULL || strlen(s) == 0) {
        return 0;
    }
    n = 1;
    for (i = 0; i < strlen(s); ++i) {
        if (s[i] == ',') {
            ++n;
        }
    }
    return n;
}
char** extractNames(char *nameList, unsigned n) {
    // n is the number of names in the name list
    if (n == 0) {
        return NULL;
    }
    char **mat = (char**) de_malloc(sizeof(char*) * n);
    unsigned index = 0;
    char *tkn = NULL;
    char *copy = de_strdup(nameList);
    tkn = strtok(copy, ",");
    while (tkn != NULL) {
        mat[index] = (char*) de_malloc(sizeof(char) * (strlen(tkn) + 1));
        strcpy(mat[index++], tkn);
        tkn = strtok(NULL, ",");
    }
    free(copy);
    return mat;
}
void destroyNameList(char **names, unsigned n) {
    for (unsigned i = 0; i < n; ++i) {
        free(names[i]);
    }
    free(names);
}
void reportNames(char **names, unsigned n) {
    for (unsigned i = 0; i < n; ++i) {
        printf("name: %s\n", names[i]);
    }
}
int main(int argc, char **argv) {
    (void) (reportNames);
    char filename[kMaxStringLength];
    char *nameList = NULL; // heap allocated by parseOptions so the list length is unbounded
    int64_t excludeBlockDegreeGT = -1;
    int64_t excludeBlockDegreeLT = -1;
    bool isInclude = true; // if 0 then we are in exclude mode. 1 is include mode.
    double maxRefNFrac = -1.;
    parseOptions(argc, argv,  filename, &nameList, &isInclude, &excludeBlockDegreeGT, &excludeBlockDegreeLT, &maxRefNFrac);
    unsigned n = countNames(nameList);
    char **names = extractNames(nameList, n);
    mafFileApi_t *mfa = maf_newMfa(filename, "r");

    filterInput(mfa, names, n, isInclude, excludeBlockDegreeGT, excludeBlockDegreeLT, maxRefNFrac);

    maf_destroyMfa(mfa);
    destroyNameList(names, n);
    free(nameList);

    return EXIT_SUCCESS;
}
