cmake_minimum_required(VERSION 3.17)
project(BamStatus)

set(CMAKE_CXX_STANDARD 14)

set(INC_DIR  /Users/zhaozhan/CLionProjects/BamStatus/htslib)
set(LINK_DIR /Users/zhaozhan/CLionProjects/BamStatus/htslib/build/lib)
include_directories(${INC_DIR})
link_directories(${LINK_DIR})

aux_source_directory(/Users/zhaozhan/CLionProjects/BamStatus/htslib/ HTSLIB)

add_executable(Config ${SAMTOOLS} ${HTSLIB} ${HTSLIB_TEXT} main.cpp Duplicate.cpp Duplicate.h BamBlock.cpp BamBlock.h block_mul.cpp Buffer.cpp Buffer.h BamStatus.cpp BamStatus.h config.h Overrepresent.cpp Overrepresent.h)

add_executable(BamStatus main.cpp )

target_link_libraries(BamStatus libhts.a)