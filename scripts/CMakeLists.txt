set (CMAKE_CXX_STANDARD 17)

add_executable(treefy treefy.C)

target_include_directories(treefy PRIVATE ${ROOT_INCLUDE_DIRS})

target_link_libraries(treefy ${ROOT_LIBRARIES})
