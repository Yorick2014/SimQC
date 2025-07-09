cmake_minimum_required(VERSION 3.14)
project(QKD_Simulator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Библиотека
add_library(pulse STATIC src/pulse.cpp)

# Заголовки для библиотеки pulse
target_include_directories(pulse PUBLIC src external)

# Программа
add_executable(qkd_simulator src/main.cpp)

# Линкуем библиотеку
target_link_libraries(qkd_simulator PRIVATE pulse)

# Компиляторные флаги
# target_compile_options(qkd_simulator PRIVATE -Wall -Wextra -Wpedantic)
