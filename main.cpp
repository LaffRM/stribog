#include <iostream>
#include "Stribog.h"

uint8_t reference_message[63] =
{
   0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37,
   0x38, 0x39, 0x30, 0x31, 0x32, 0x33, 0x34, 0x35,
   0x36, 0x37, 0x38, 0x39, 0x30, 0x31, 0x32, 0x33,
   0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x30, 0x31,
   0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39,
   0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37,
   0x38, 0x39, 0x30, 0x31, 0x32, 0x33, 0x34, 0x35,
   0x36, 0x37, 0x38, 0x39, 0x30, 0x31, 0x32
};

uint8_t reference_message_hash[64] =
{
	0x48, 0x6f, 0x64, 0xc1, 0x91, 0x78, 0x79, 0x41, 0x7f, 0xef, 0x08, 0x2b, 0x33, 0x81, 0xa4, 0xe2,
	0x11, 0xc3, 0x24, 0xf0, 0x74, 0x65, 0x4c, 0x38, 0x82, 0x3a, 0x7b, 0x76, 0xf8, 0x30, 0xad, 0x00,
	0xfa, 0x1f, 0xba, 0xe4, 0x2b, 0x12, 0x85, 0xc0, 0x35, 0x2f, 0x22, 0x75, 0x24, 0xbc, 0x9a, 0xb1,
	0x62, 0x54, 0x28, 0x8d, 0xd6, 0x86, 0x3d, 0xcc, 0xd5, 0xb9, 0xf5, 0x4a, 0x1a, 0xd0, 0x54, 0x1b
};


int main()
{
	char p;
	uint8_t hash_t[STRIBOG_HASH_LENGTH] = {'\0'};
	uint8_t hash_test[STRIBOG_HASH_LENGTH] = { '\0' };
	Stribog a;

	a.stribog_hash(hash_t, reference_message, sizeof(reference_message));

	a.hash_reverce(hash_t);

	std::cout << hash_t << std::hex << std::endl;

	if (!memcmp(hash_t, reference_message_hash, STRIBOG_HASH_LENGTH))
	{
		std::cout <<"hash Successful"<<std::endl;
	}

	else
	{
		std::cout << "hash Failed" << std::endl;
	}

	return 0;
}