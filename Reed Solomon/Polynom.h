#ifndef POLYNOM_H
#define POLYNOM_H

#include <stdint.h>
#include <string.h>
#include <vector>

using namespace std;

namespace RS {

	struct Poly {
		Poly() {
			coefficents_ = std::vector<uint8_t>(0);
		}
		

		/* @brief Fügt Nummer am Ende des Polynoms hinzu */
		void Append(uint8_t num) {
			coefficents_.push_back(num);
		}

		/* @brief Initialisieren des Polynoms mit Startwerten auf 0*/
		void Init(uint8_t length) {
			this->coefficents_ = std::vector<uint8_t>(length);
			Reset();
		}

		/* @brief Polynom wird auf Null gesetzt*/
		void Reset() {
			for (int i = 0; i < coefficents_.size(); i++) {
				coefficents_.at(i) = 0;
			}
		}

		void Set(uint8_t* array, uint8_t count) {
			for (uint8_t i = 0; i < count; i++) {
				coefficents_.at(i) = array[i];
			}
		}

		#define poly_max(a, b) ((a > b) ? (a) : (b))

		void Copy(Poly tmp, int end, int begin = 0) {
			Init(0);
			for (int i = begin; i < end; i++) {
				Append(tmp.at(i));
			}
		}

		uint8_t& at(uint8_t i) {
			return coefficents_.at(i);
		}

		uint8_t length() {
			return coefficents_.size();
		}

	protected:
		std::vector<uint8_t> coefficents_;
	};
};

#endif 