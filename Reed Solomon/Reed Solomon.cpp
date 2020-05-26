/*  Diese Implementation basiert auf der Implementation von https://github.com/mersinvald/Reed-Solomon,
    sowie dem Beispiel-Code von https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders#RS_encoding
*/

#include <iostream>
#include <fstream>
#include <exception>
#include <stdint.h>
#include <time.h>
#include <ctime>
#include <chrono>

#include "RS_Codec.h"

constexpr int MESSAGE_LENGTH = 24;
constexpr int ECC_LENGTH = 12;

std::vector<std::string> readFile(std::string filePath) {
    std::ifstream file(filePath);
    std::vector<std::string> data;

    if (file.is_open()) {

        while (!file.eof()) {

            char inp[(MESSAGE_LENGTH + 1)] = "";

            file.read(inp, MESSAGE_LENGTH);

            std::string tmp = inp;

            if (tmp.length() < MESSAGE_LENGTH) {
                for (int i = tmp.length(); i < MESSAGE_LENGTH; i++) {
                    tmp.push_back(0);
                }
            }
            data.push_back(tmp);
        }
    }
    else {
        throw std::runtime_error("Falscher Pfad");
    }

    std::cout << "Daten eingelesen!" << std::endl;

    return data;
}

int main()
{
    srand(time(NULL));

    try {
        RS::RS_Codec<MESSAGE_LENGTH, ECC_LENGTH> RS;

        std::vector<std::string> data;// = readFile("C:\\Users\\Domin\\Desktop\\test.txt");

        data.push_back("Eine wichtige Nachricht.");

        auto start = std::chrono::high_resolution_clock::now();

        RS.encode(data);

        std::cout << data.at(0) << std::endl;

        auto elapsed = std::chrono::high_resolution_clock::now() - start;

        RS.corrupt(data, 7);

        std::cout << data.at(0) << std::endl;

        auto start2 = std::chrono::high_resolution_clock::now();

        uint8_t array[] = {0, 1, 2};

        RS.decode(data, array, 3);

        std::cout << data.at(0) << std::endl;

        auto elapsed2 = std::chrono::high_resolution_clock::now() - start2;

        long long encode_time = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
        long long decode_time = std::chrono::duration_cast<std::chrono::microseconds>(elapsed2).count();

        std::cout << "Anzahl Bloecke: " << data.size() << std::endl;
        std::cout << "Encode-Zeit:  " << encode_time << " mykro_sekunden" << std::endl;
        std::cout << "Decode-Zeit:  " << decode_time << " mykro_sekunden" << std::endl;
    }
    catch (runtime_error e) {
        std::cout << e.what() << std::endl;
    }
}