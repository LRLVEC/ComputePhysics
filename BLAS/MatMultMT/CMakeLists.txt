project(MatMultMT)

# set(CMAKE_THREAD_PREFER_PTHREAD ON)
# set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads)
add_executable(MatMultMT MatMultMT.cpp)
target_link_libraries(MatMultMT Threads::Threads)