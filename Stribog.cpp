#include "Stribog.h"



Stribog::Stribog()
{
	s_data = new stribog_data;
}


Stribog::~Stribog()
{
}

void Stribog::xor_vect(const uint8_t *a, const uint8_t *b, uint8_t *c)
{
	int i = 0;

	for (i = 0; i < BLOCK_SIZE; i++)
	{
		c[i] = a[i] ^ b[i];
	}
}

void Stribog::add_512(const uint8_t *a, const uint8_t *b, uint8_t *c)
{
	uint16_t tmp = 0;
	for (int i = BLOCK_SIZE - 1; i >= 0; i--)
	{
		tmp = (uint16_t)a[i] + (uint16_t)b[i] + (tmp >> 8);
		c[i] = (uint8_t)(tmp & 0xff);
	}
}

void Stribog::sub_s(uint8_t *vect)
{
	uint8_t tmp[BLOCK_SIZE];

	for (int i = BLOCK_SIZE - 1; i >= 0; i--)
	{
		tmp[i] = sbox_pi[vect[i]];
	}

	memcpy(vect, tmp, BLOCK_SIZE);
}

void Stribog::perm_p(uint8_t *vect)
{
	uint8_t tmp[BLOCK_SIZE] = {'0'};

	for (int i = BLOCK_SIZE - 1; i >= 0; i--)
	{
		tmp[i] = vect[tau[i]];
	}

	memcpy(vect, tmp, BLOCK_SIZE);
}

void Stribog::linr_l(uint8_t *vect)
{
	uint64_t* tmp_in = (uint64_t *)vect;
	uint64_t tmp_out[8];

	memset(tmp_out, 0x00, BLOCK_SIZE);

	for (int i = 7; i >= 0; i--)
	{
		for (int j = 63; j >= 0; j--)
		{
			if ((tmp_in[i] >> j) & 1)
			{
				tmp_out[i] = tmp_out[i] ^ a_mtrx[63 - j];
			}
		}
	}

	memcpy(vect, tmp_out, BLOCK_SIZE);

}

void Stribog::spl(uint8_t *vect)
{
	sub_s(vect);
	perm_p(vect);
	linr_l(vect);
}
void Stribog::get_key(uint8_t *K, int i)
{
	xor_vect(K, C[i], K);
	spl(K);
}

void Stribog::e(uint8_t *K, const uint8_t *m, uint8_t *vect)
{
	memcpy(K, K, BLOCK_SIZE);
	xor_vect(m, K, vect);

	for (int i = 0; i < 12; i++)
	{
		spl(vect);
		get_key(K, i);
		xor_vect(vect, K, vect);
	}
}

void Stribog::g(uint8_t *h, uint8_t *N, const uint8_t *m)
{
	uint8_t K[BLOCK_SIZE] = {'\0'};
	uint8_t tmp[BLOCK_SIZE] = {'\0'};

	xor_vect(N, h, K);

	spl(K);

	e(K, m, tmp);

	xor_vect(tmp, h, tmp);
	xor_vect(tmp, m, h);
}

void Stribog::padding(stribog_data *sdt)
{
	uint8_t tmp[BLOCK_SIZE];
	if (sdt->buf_size < BLOCK_SIZE)
	{
		memset(tmp, 0x00, BLOCK_SIZE);
		memcpy(tmp, sdt->buffer, sdt->buf_size);

		tmp[sdt->buf_size] = 0x01;
		memcpy(sdt->buffer, tmp, BLOCK_SIZE);
	}
}

void Stribog::init(stribog_data *sdt)
{
	memset(sdt, 0x00, sizeof(stribog_data));
	memcpy(sdt->h, iv512, BLOCK_SIZE);
	sdt->v_512[1] = 0x02;
}

void Stribog::step_2(stribog_data *sdt, uint8_t *data)
{

	g(sdt->h, sdt->N, data);
	add_512(sdt->N, sdt->v_512, sdt->N);
	add_512(sdt->S, data, sdt->S);

}

void Stribog::step_3(stribog_data *sdt)
{
	uint8_t tmp[BLOCK_SIZE];
	memset(tmp, 0x00, BLOCK_SIZE);

	tmp[1] = ((sdt->buf_size * 8) >> 8) & 0xff;
	tmp[0] = (sdt->buf_size * 8) & 0xff;

	padding(sdt);

	g(sdt->h, sdt->N, (const uint8_t*)&(sdt->buffer));

	add_512(sdt->N, tmp, sdt->N);
	add_512(sdt->S, sdt->buffer, sdt->S);

	g(sdt->h, sdt->v_0, (uint8_t*)&(sdt->N));
	g(sdt->h, sdt->v_0, (uint8_t*)&(sdt->S));

	memcpy(sdt->hash, sdt->h, BLOCK_SIZE);
}

void Stribog::hash_update(stribog_data *sdt, uint8_t *data, size_t data_len)
{
	size_t empt_size;

	while ((data_len > 63) && (sdt->buf_size) == 0)
	{
		step_2(sdt, data);
		data = data + 64;
		data_len = data_len - 64;
	}

	while (data_len)
	{
		empt_size = 64 - sdt->buf_size;
		if (empt_size > data_len)
		{
			empt_size = data_len;
		}

		memcpy(&sdt->buffer[sdt->buf_size], data, empt_size);
		sdt->buf_size = sdt->buf_size + empt_size;
		data_len = data_len - empt_size;
		data = data + empt_size;
		if (sdt->buf_size == 64)
		{
			step_2(sdt, sdt->buffer);
			sdt->buf_size = 0;
		}
	}
}

void Stribog::get_hash(stribog_data *sdt)
{
	step_3(sdt);
	sdt->buf_size = 0;
}
void Stribog::stribog_hash(uint8_t *hash, const uint8_t *data, size_t data_len)
{
	init(s_data);

	uint8_t *data_tmp = new uint8_t [data_len];

	memcpy(data_tmp, data, data_len);

	hash_update(s_data, data_tmp, data_len);

	get_hash(s_data);

	memcpy(hash, s_data->hash, BLOCK_SIZE);
}

void Stribog::hash_reverce(uint8_t *hash)
{
	uint8_t hash_reverce[STRIBOG_HASH_LENGTH];

	for (int i = 0; i < STRIBOG_HASH_LENGTH; i++)
	{
		hash_reverce[i] = hash[STRIBOG_HASH_LENGTH - 1 - i];
	}

	memcpy(hash, hash_reverce, STRIBOG_HASH_LENGTH);
}