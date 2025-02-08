#include <gtest/gtest.h>
#include <fstream>
#include <string>

// This test verifies that the unit test framework is available and
// that we can access files in the "TestFiles" directory.
TEST(ExampleTest, ExampleTest) {
    const std::string filePath = "TestFiles/ExampleFile.txt"; 
    std::ifstream file(filePath);

    // Verify the file was successfully opened
    ASSERT_TRUE(file.is_open()) << "Failed to open file: " << filePath;

    // Read the contents of the file
    std::string content;
    std::getline(file, content);

    // Verify that the file is not empty (example assumption that the file should contain some text)
    EXPECT_TRUE(content.find("Example text") != std::string::npos);
    
    // Close the file
    file.close();
}
