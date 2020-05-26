#ifndef RS_CODEC_H
#define RS_CODEC_H


#include <string>	
#include <iostream>
#include <exception>
#include "Galois_Field.h"

namespace RS {
	template <const int msg_length,
			  const int ecc_length>
		class RS_Codec
	{
		int block_length;
		Poly Generator;

	public:
		RS_Codec() {
			this->block_length = msg_length + ecc_length;
			Generator = GeneratorPoly();
		}

		void encode(std::vector<std::string>& data) {
			// Überprüfung der Länge
			if (data[0].size() > msg_length && block_length > 255) {
				throw std::runtime_error("Blockgröße ist zu groß.");
			}

			// Jeder einzelne Block wird kodiert
			for (std::string& string : data) {
				encode_Blocks(string);
			}
		}

		void corrupt(std::vector<std::string>& data, unsigned int amount) {
			for (std::string& current : data) {
				for (int i = 0; i < amount; i++) {
					current.at(i) = 0;
				}
			}
		}

		void decode(std::vector<std::string>& data, uint8_t* erase_pos = NULL, size_t erase_count = 0) {
			if (data[0].size() > msg_length && block_length > 255) {
				throw std::runtime_error("Blockgröße ist zu groß.");
			}

			for (int i = 0; i < data.size(); i++) {
				decode_Blocks(data.at(i), erase_pos, erase_count);
			}
		}

	private:

		void dataToPoly(std::string& data, Poly& msg) {
			for (int c = 0; c < block_length; c++) {
				msg.Append((uint8_t)data.c_str()[c]);
			}
		}

		void polyToData(std::string& data, Poly& msg) {
			std::string tmp;
			for (int i = 0; i < msg_length; i++) {
				tmp.push_back(msg.at(i));
			}

			data = tmp;
		}

		void encode_Blocks(std::string& data) {

			// Nachricht wird in ein Polynom gepackt
			Poly msg;
			for (int c = 0; c < block_length; c++) {
				if (c < data.size()) {
					msg.Append((uint8_t)data.c_str()[c]);
				}
				else {
					msg.Append(0);
				}
			}

			// Bestimmen der Redundanz
			uint8_t coef = 0;
			for (uint8_t i = 0; i < msg_length; i++) {
				coef = msg.at(i);
				if (coef != 0) {
					for (uint32_t j = 1; j < Generator.length(); j++) {
						msg.at(i + j) ^= GF::mul(Generator.at(j), coef);
					}
				}
			}

			std::cout << "Rest:" << std::endl;
			// Redundanz an das Ende der Nachricht anfuegen
			std::string tmp;
			for (int c = msg_length; c < block_length; c++) {
				std::cout << (int)msg.at(c) << "|";
				tmp.push_back(msg.at(c));
			}
			std::cout << std::endl;
			data += tmp;
		}

		void decode_Blocks(std::string& data, uint8_t* erase_pos = NULL, size_t erase_count = 0) {

			// Daten in Polynom schreiben
			Poly msg;
			dataToPoly(data, msg);

			Poly epos;
			if (erase_pos == NULL) {
				epos.Init(0);
			}
			else {
				epos.Init(erase_count);
				epos.Set(erase_pos, erase_count);

				for (uint8_t i = 0; i < epos.length(); i++) {
					msg.at(epos.at(i)) = 0;
				}
			}

			// Zuviele Fehler
			if (epos.length() > ecc_length) {
				throw runtime_error("Zuviele Fehler!");
			}

			Poly synd, eloc, reloc, err, forney;

			// Syndrome berechnen
			synd = calc_syndrom(msg);

			// Nach Fehler schauen
			bool has_errors = false;
			for (uint8_t i = 0; i < synd.length(); i++) {
				if (synd.at(i) != 0) {
					has_errors = true;
					break;
				}
			}

			// Falls keine Fehler vorhanden
			if (!has_errors) {
				polyToData(data, msg);
				return;
			}

			forney = calc_forney_syndromes(synd, epos);
			eloc = find_error_locator(forney, NULL, epos.length()); 

			// Fehlerpolynom-Ausgabe
			std::cout << "Fehlerpolynom" << std::endl;
			for (int i = 0; i < eloc.length(); i++) {
				std::cout << (int)eloc.at(i) << "|";
			}
			std::cout << std::endl;

			// Syndrom umkehren
			reloc.Init(eloc.length());
			for (int8_t i = eloc.length() - 1, j = 0; i >= 0; i--, j++) {
				reloc.at(j) = eloc.at(i);
			}

			// Fehler finden
			bool ok;
			err = find_errors(reloc, ok);
			if (!ok || (err.length() == 0 && erase_count == 0)) { 
				throw runtime_error("Es kam zu einem Fehler bei der Fehlersuche!");
			}

			// Fehlerpolynom-Ausgabe
			std::cout << "Fehlerstellen" << std::endl;
			for (int i = 0; i < err.length(); i++) {
				std::cout << (int)err.at(i) << "|";
			}
			std::cout << std::endl;

			// Hinzufuegen der gefundenen Fehler
			for (uint8_t i = 0; i < err.length(); i++) {
				epos.Append(err.at(i));
			}

			// Korrigieren der Fehler
			msg = correct_errata(synd, epos, msg);

			polyToData(data, msg);
		}

		Poly GeneratorPoly() {
			Poly gen;
			gen.Init(1);
			gen.at(0) = 1;

			Poly mul;
			mul.Init(2);

			for (int i = 0; i < ecc_length; i++) {
				mul.at(0) = 1;
				mul.at(1) = GF::pow(2, i);
				gen = GF::poly_mul(gen, mul);
			}

			return gen;
		}

		Poly calc_syndrom(Poly& msg) {
			Poly synd;
			synd.Init(ecc_length + 1);
			synd.at(0) = 0;

			for (uint8_t i = 1; i < ecc_length + 1; i++) {
				synd.at(i) = GF::poly_eval(msg, GF::pow(2, i - 1));
			}

			return synd;
		}

		Poly calc_forney_syndromes(Poly& synd, Poly& erasures_pos) {
			Poly erase_pos_reversed;
			Poly fsynd;

			for (uint8_t i = 0; i < erasures_pos.length(); i++) {
				erase_pos_reversed.Append(block_length - 1 - erasures_pos.at(i));
			}

			fsynd.Init(synd.length() - 1);
			fsynd.Copy(synd, synd.length(), 1);

			uint8_t x;
			for (uint8_t i = 0; i < erasures_pos.length(); i++) {
				x = GF::pow(2, erase_pos_reversed.at(i));
				for (uint8_t j = 0; j < fsynd.length() - 1; j++) {
					fsynd.at(j) = GF::mul(fsynd.at(j), x) ^ fsynd.at(j + 1);
				}
			}

			return fsynd;
		}
		
		Poly find_error_locator(Poly& synd, Poly* erase_loc = NULL, uint8_t erase_count = 0) {
			Poly error_loc;
			Poly err_loc;
			Poly old_loc;
			Poly temp;
			Poly temp2;
			
			if (erase_loc != NULL) {
				err_loc.Copy(*erase_loc, erase_loc->length());
				old_loc.Copy(*erase_loc, erase_loc->length());
			}
			else {
				err_loc.Init(1);
				old_loc.Init(1);
				err_loc.at(0) = 1;
				old_loc.at(0) = 1;
			}

			uint8_t synd_shift = 0;
			if (synd.length() > ecc_length) {
				synd_shift = synd.length() - ecc_length;
			}

			uint8_t K = 0;
			uint8_t delta = 0;
			uint8_t index;

			for (uint8_t i = 0; i < ecc_length - erase_count; i++) {
				if (erase_loc != NULL) {
					K = erase_count + i + synd_shift;
				}
				else {
					K = i + synd_shift;
				}

				delta = synd.at(K);
				
				for (uint8_t j = 1; j < err_loc.length(); j++) {
					index = err_loc.length() - j - 1;
					delta ^= GF::mul(err_loc.at(index), synd.at(K - j));
				}

				old_loc.Append(0);

				if (delta != 0) {
					if (old_loc.length() > err_loc.length()) {
						temp = GF::poly_scale(old_loc, delta);
						old_loc = GF::poly_scale(err_loc, GF::inverse(delta));
						err_loc.Copy(temp, temp.length());
					}
					temp = GF::poly_scale(old_loc, delta);
					temp2 = GF::poly_add(err_loc, temp);
					err_loc.Copy(temp2, temp2.length());
				}
			}

			uint32_t shift = 0;
			while (err_loc.length() && err_loc.at(shift) == 0) {
				shift++;
			}

			uint32_t errs = err_loc.length() - 1;

			if ((int)((errs - erase_count) * 2 + erase_count) > ecc_length) {
				throw runtime_error("Nachricht kann aufgrund zu vieler Fehler nicht korrigiert werden!");
			}

			err_loc.Copy(err_loc, err_loc.length(), shift);
			return err_loc;
		}

		Poly find_errors(Poly& err_loc, bool& ok) {
			Poly err_pos;
			err_pos.Init(0);

			uint8_t errs = err_loc.length() - 1;
			for (uint8_t i = 0; i < block_length; i++) {
				if (GF::poly_eval(err_loc, GF::pow(2, i)) == 0) {
					err_pos.Append(block_length - 1 - i);
				}
			}

			if (err_pos.length() != errs) {
				throw runtime_error("Konnte Fehlerstellen nicht finden.");
			}

			ok = true;
			return err_pos;
		}
		
		Poly correct_errata(Poly& synd, Poly& err_pos, Poly& msg) {
			Poly coef_pos;
			coef_pos.Init(err_pos.length());

			for (uint8_t i = 0; i < err_pos.length(); i++) {
				coef_pos.at(i) = msg.length() - 1 - err_pos.at(i);
			}

			Poly errata_loc;
			errata_loc = find_errata_locator(coef_pos);

			// Umgekehrtes Syndrom-Polynom
			Poly rsynd;
			rsynd.Init(synd.length());

			for (int8_t i = synd.length() - 1, j = 0; i >= 0; i--, j++) {
				rsynd.at(j) = synd.at(i);
			}

			// Umgedrehtes eval Polynom
			Poly re_eval;
			re_eval = find_error_evaluator(rsynd, errata_loc, errata_loc.length() - 1);

			// Zurueck drehen
			Poly e_eval;
			e_eval.Init(re_eval.length());
			for (int8_t i = re_eval.length() - 1, j = 0; i >= 0; i--, j++) {
				e_eval.at(j) = re_eval.at(i);
			}

			// Speichert die Fehlerpositionen
			Poly X;
			X.Init(0);

			uint16_t l;
			for (uint8_t i = 0; i < coef_pos.length(); i++) {
				l = 255 - coef_pos.at(i);
				X.Append(GF::pow(2, -l));
			}

			// Magnitude Polynom
			Poly E;
			E.Init(msg.length());

			uint8_t Xi_inv;

			Poly err_loc_prime_temp;

			uint8_t err_loc_prime;
			uint8_t y;

			for (uint8_t i = 0; i < X.length(); i++) {
				Xi_inv = GF::inverse(X.at(i));

				err_loc_prime_temp.Init(0);
				for (uint8_t j = 0; j < X.length(); j++) {
					if (j != i) {
						err_loc_prime_temp.Append(GF::sub(1, GF::mul(Xi_inv, X.at(j))));
					}
				}

				err_loc_prime = 1;
				for (uint8_t j = 0; j < err_loc_prime_temp.length(); j++) {
					err_loc_prime = GF::mul(err_loc_prime, err_loc_prime_temp.at(j));
				}

				y = GF::poly_eval(re_eval, Xi_inv);
				y = GF::mul(GF::pow(X.at(i), 1), y);

				E.at(err_pos.at(i)) = GF::div(y, err_loc_prime);
			}

			msg = GF::poly_add(msg, E);

			return msg;
		}
		
		Poly find_errata_locator(Poly& e_pos) {
			Poly e_loc;
			e_loc.Init(1);
			e_loc.at(0) = 1;

			Poly add1;
			add1.Init(1);

			Poly add2;
			add2.Init(2);

			for (int i = 0; i < e_pos.length(); i++) {
				add1.at(0) = 1;
				add2.at(0) = GF::pow(2, e_pos.at(i));
				add2.at(1) = 0;
				e_loc = GF::poly_mul(e_loc, GF::poly_add(add1, add2));
			}
			return e_loc;
		}

		Poly find_error_evaluator(Poly& synd, Poly& errata_loc, uint8_t ecclen) {
			Poly remainder;
			Poly mulp;
			Poly divisor;
			divisor.Init(ecclen + 2);
			divisor.at(0) = 1;

			mulp = GF::poly_mul(synd, errata_loc);

			remainder = GF::poly_div(mulp, divisor);

			return remainder;
		}
	};
};
#endif 