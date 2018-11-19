#ifndef R3BCALIFAFEBEXREADER_H
#define R3BCALIFAFEBEXREADER_H

#include "R3BReader.h"
class TClonesArray;
class FairLogger;

struct EXT_STR_h101_CALIFA_t;
typedef struct EXT_STR_h101_CALIFA_t EXT_STR_h101_CALIFA;
class ext_data_struct_info;
/**
 * A reader of CALIFA FEBEX data with UCESB.
 * Receives mapped raw data and converts it to R3BRoot objects.
 * @author H. Alvarez
 * @since May 16, 2017
 */
class R3BCalifaFebexReader : public R3BReader {
	public:
		/**
		 * Standard constructor.
		 * Creates instance of the reader. To be called in the steering macro.
		 * @param Pointer to a full C structure generated by the Ucesb unpacker.
		 */
		R3BCalifaFebexReader(EXT_STR_h101_CALIFA *, UInt_t);
		~R3BCalifaFebexReader();

		Bool_t Init(ext_data_struct_info *);
		Bool_t Read();
		void Reset();

	private:
		/* An event counter */
		unsigned int fNEvent;
		/* Reader specific data structure from ucesb */
		EXT_STR_h101_CALIFA* fData;
		/* Data offset */
		UInt_t fOffset;
		/* FairLogger */
		FairLogger*	fLogger;

                Bool_t fOnline;//Don't store data for online

	  TClonesArray* fArray; /**< Output array. */

	public:
		ClassDef(R3BCalifaFebexReader, 0);
};

#endif
