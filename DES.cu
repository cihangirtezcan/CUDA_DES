#include <Windows.h>
#include <stdio.h>
#include <cuda.h>
#include <time.h>
#include "tables_PRESENT.inc"
#include "tables_DES.inc"
#define Exhaustive 65536
#define THREAD 1024
#define BLOCK 512

bit64 ciphertext[THREAD * BLOCK] = { 0 };


__global__ void DES_CTR(bit64* subkey_d, bit64* ciphertext_d, bit64* expansion_table0_d, bit32* expansion_table1_d, bit64* expansion_table2_d, bit64* expansion_table3_d, bit32* s_permutation_table0_d, bit32* s_permutation_table1_d, bit32* s_permutation_table2_d, bit32* s_permutation_table3_d, bit32* s_permutation_table4_d, bit32* s_permutation_table5_d, bit32* s_permutation_table6_d, bit32* s_permutation_table7_d) {
	__shared__ bit64 expansion_table0[256], expansion_table2[256], expansion_table3[256], subkey[16];
	__shared__ bit32 expansion_table1[256], s_permutation_table0[64], s_permutation_table1[64], s_permutation_table2[64], s_permutation_table3[64], s_permutation_table4[64], s_permutation_table5[64], s_permutation_table6[64], s_permutation_table7[64];
	if (threadIdx.x < 256) {
		expansion_table0[threadIdx.x] = expansion_table0_d[threadIdx.x];
		expansion_table1[threadIdx.x] = expansion_table1_d[threadIdx.x];
		expansion_table2[threadIdx.x] = expansion_table2_d[threadIdx.x];
		expansion_table3[threadIdx.x] = expansion_table3_d[threadIdx.x];
	}
	if (threadIdx.x < 64) {
		s_permutation_table0[threadIdx.x] = s_permutation_table0_d[threadIdx.x];
		s_permutation_table1[threadIdx.x] = s_permutation_table1_d[threadIdx.x];
		s_permutation_table2[threadIdx.x] = s_permutation_table2_d[threadIdx.x];
		s_permutation_table3[threadIdx.x] = s_permutation_table3_d[threadIdx.x];
		s_permutation_table4[threadIdx.x] = s_permutation_table4_d[threadIdx.x];
		s_permutation_table5[threadIdx.x] = s_permutation_table5_d[threadIdx.x];
		s_permutation_table6[threadIdx.x] = s_permutation_table6_d[threadIdx.x];
		s_permutation_table7[threadIdx.x] = s_permutation_table7_d[threadIdx.x];
	}
	if (threadIdx.x < 16) subkey[threadIdx.x] = subkey_d[threadIdx.x];
	__syncthreads();
	bit64 temp_exp, temp;
	bit32 plaintext_left, plaintext_right;
	bit64 threadIndex = blockIdx.x * blockDim.x + threadIdx.x;
	int i, j;
	for (i = 0; i < Exhaustive; i++) {
		plaintext_left = 0;
		plaintext_right = threadIndex + i * THREAD * BLOCK;
		for (j = 0; j < 16; j++) {
			temp_exp = expansion_table1[(plaintext_right >> 8) & 0xff] ^ expansion_table0[plaintext_right & 0xff] ^ expansion_table2[(plaintext_right >> 16) & 0xff] ^ expansion_table3[plaintext_right >> 24] ^ subkey[j];
			plaintext_left ^= s_permutation_table0[temp_exp & 0x3f] ^ s_permutation_table1[(temp_exp >> 6) & 0x3f] ^ s_permutation_table2[(temp_exp >> 12) & 0x3f] ^ s_permutation_table3[(temp_exp >> 18) & 0x3f] ^ s_permutation_table4[(temp_exp >> 24) & 0x3f] ^ s_permutation_table5[(temp_exp >> 30) & 0x3f] ^ s_permutation_table6[(temp_exp >> 36) & 0x3f] ^ s_permutation_table7[temp_exp >> 42];
			temp = plaintext_left;
			plaintext_left = plaintext_right;
			plaintext_right = temp;
		}
		temp = plaintext_left; temp = temp << 32; temp ^= plaintext_right;
		ciphertext_d[threadIndex] = temp;
	}
}
__global__ void DES_exhaustive(bit64 *plaintext_d, bit64 key2, bit32 plaintext_left2,bit32 plaintext_right2,bit32 ciphertext_left2,bit32 ciphertext_right2,bit64 *expansion_table0_d,bit32 *expansion_table1_d,bit64 *expansion_table2_d,bit64 *expansion_table3_d,bit32 *PC2_table0_d,bit32 *PC2_table1_d,bit32 *PC2_table2_d,bit64 *PC2_table3_d,bit64 *PC2_table4_d,bit64 *PC2_table5_d,bit64 *PC2_table6_d,bit32 *s_permutation_table0_d,bit32 *s_permutation_table1_d,bit32 *s_permutation_table2_d,bit32 *s_permutation_table3_d,bit32 *s_permutation_table4_d,bit32 *s_permutation_table5_d,bit32 *s_permutation_table6_d,bit32 *s_permutation_table7_d){
	__shared__ bit64 expansion_table0[256],expansion_table2[256],expansion_table3[256];
	__shared__ bit64 PC2_table3[256],PC2_table4[256],PC2_table5[256],PC2_table6[256];
	__shared__ bit32 expansion_table1[256],PC2_table0[256],PC2_table1[256],PC2_table2[256],s_permutation_table0[64],s_permutation_table1[64],s_permutation_table2[64],s_permutation_table3[64],s_permutation_table4[64],s_permutation_table5[64],s_permutation_table6[64],s_permutation_table7[64];
	if (threadIdx.x < 256) {
		expansion_table0[threadIdx.x] = expansion_table0_d[threadIdx.x];
		expansion_table1[threadIdx.x] = expansion_table1_d[threadIdx.x];
		expansion_table2[threadIdx.x] = expansion_table2_d[threadIdx.x];
		expansion_table3[threadIdx.x] = expansion_table3_d[threadIdx.x];
		PC2_table0[threadIdx.x] = PC2_table0_d[threadIdx.x];
		PC2_table1[threadIdx.x] = PC2_table1_d[threadIdx.x];
		PC2_table2[threadIdx.x] = PC2_table2_d[threadIdx.x];
		PC2_table3[threadIdx.x] = PC2_table3_d[threadIdx.x];
		PC2_table4[threadIdx.x] = PC2_table4_d[threadIdx.x];
		PC2_table5[threadIdx.x] = PC2_table5_d[threadIdx.x];
		PC2_table6[threadIdx.x] = PC2_table6_d[threadIdx.x];
	}
	if (threadIdx.x < 64) {
		s_permutation_table0[threadIdx.x] = s_permutation_table0_d[threadIdx.x];
		s_permutation_table1[threadIdx.x] = s_permutation_table1_d[threadIdx.x];
		s_permutation_table2[threadIdx.x] = s_permutation_table2_d[threadIdx.x];
		s_permutation_table3[threadIdx.x] = s_permutation_table3_d[threadIdx.x];
		s_permutation_table4[threadIdx.x] = s_permutation_table4_d[threadIdx.x];
		s_permutation_table5[threadIdx.x] = s_permutation_table5_d[threadIdx.x];
		s_permutation_table6[threadIdx.x] = s_permutation_table6_d[threadIdx.x];
		s_permutation_table7[threadIdx.x] = s_permutation_table7_d[threadIdx.x];
	}
	__syncthreads();
	bit32 plaintext_leftr=plaintext_left2, plaintext_rightr=plaintext_right2, ciphertext_left=ciphertext_left2, ciphertext_right=ciphertext_right2;
	bit64 key_real=key2+blockIdx.x*blockDim.x+threadIdx.x,temp_exp,subkey,key;
	bit32 plaintext_left,plaintext_right;
	int i,j;
	key=key_real;
	for (i=0;i<Exhaustive;i++) { 
		plaintext_left=plaintext_leftr;
		plaintext_right=plaintext_rightr;
		for (j=1;j<9;j++) {
			if (j==1 || j==5) key=((key<<1)&0xffffffeffffffe)^((key>>27)&0x10000001);
			else key=((key<<2)&0xffffffcffffffc)^((key>>26)&0x30000003);
			subkey=PC2_table0[key&0xff]^PC2_table1[(key>>8)&0xff]^PC2_table2[(key>>16)&0xff]^PC2_table3[(key>>24)&0xff]^PC2_table4[(key>>32)&0xff]^PC2_table5[(key>>40)&0xff]^PC2_table6[key>>48];
			temp_exp=expansion_table1[(plaintext_right>>8)&0xff]^expansion_table0[plaintext_right&0xff]^expansion_table2[(plaintext_right>>16)&0xff]^expansion_table3[plaintext_right>>24]^subkey;
			plaintext_left^=s_permutation_table0[temp_exp&0x3f]^s_permutation_table1[(temp_exp>>6)&0x3f]^s_permutation_table2[(temp_exp>>12)&0x3f]^s_permutation_table3[(temp_exp>>18)&0x3f]^s_permutation_table4[(temp_exp>>24)&0x3f]^s_permutation_table5[(temp_exp>>30)&0x3f]^s_permutation_table6[(temp_exp>>36)&0x3f]^s_permutation_table7[temp_exp>>42];
			if (j==1) {
				key=((key<<1)&0xffffffeffffffe)^((key>>27)&0x10000001);
				subkey=PC2_table0[key&0xff]^PC2_table1[(key>>8)&0xff]^PC2_table2[(key>>16)&0xff]^PC2_table3[(key>>24)&0xff]^PC2_table4[(key>>32)&0xff]^PC2_table5[(key>>40)&0xff]^PC2_table6[key>>48];
				temp_exp=expansion_table0[plaintext_left&0xff]^expansion_table1[(plaintext_left>>8)&0xff]^expansion_table2[(plaintext_left>>16)&0xff]^expansion_table3[plaintext_left>>24]^subkey;
				plaintext_right^=s_permutation_table0[temp_exp&0x3f]^s_permutation_table1[(temp_exp>>6)&0x3f]^s_permutation_table2[(temp_exp>>12)&0x3f]^s_permutation_table3[(temp_exp>>18)&0x3f]^s_permutation_table4[(temp_exp>>24)&0x3f]^s_permutation_table5[(temp_exp>>30)&0x3f]^s_permutation_table6[(temp_exp>>36)&0x3f]^s_permutation_table7[temp_exp>>42];

			}
			else if (j==8) {
				if (plaintext_left==ciphertext_right) {
					key=((key<<1)&0xffffffeffffffe)^((key>>27)&0x10000001);
					subkey=PC2_table0[key&0xff]^PC2_table1[(key>>8)&0xff]^PC2_table2[(key>>16)&0xff]^PC2_table3[(key>>24)&0xff]^PC2_table4[(key>>32)&0xff]^PC2_table5[(key>>40)&0xff]^PC2_table6[key>>48];
					temp_exp=expansion_table1[(plaintext_left>>8)&0xff]^expansion_table0[plaintext_left&0xff]^expansion_table2[(plaintext_left>>16)&0xff]^expansion_table3[plaintext_left>>24]^subkey;
					plaintext_right^=s_permutation_table0[temp_exp&0x3f]^s_permutation_table1[(temp_exp>>6)&0x3f]^s_permutation_table2[(temp_exp>>12)&0x3f]^s_permutation_table3[(temp_exp>>18)&0x3f]^s_permutation_table4[(temp_exp>>24)&0x3f]^s_permutation_table5[(temp_exp>>30)&0x3f]^s_permutation_table6[(temp_exp>>36)&0x3f]^s_permutation_table7[temp_exp>>42];
				}
			}
			else {
				key=((key<<2)&0xffffffcffffffc)^((key>>26)&0x30000003);
				subkey=PC2_table0[key&0xff]^PC2_table1[(key>>8)&0xff]^PC2_table2[(key>>16)&0xff]^PC2_table3[(key>>24)&0xff]^PC2_table4[(key>>32)&0xff]^PC2_table5[(key>>40)&0xff]^PC2_table6[key>>48];
				temp_exp=expansion_table1[(plaintext_left>>8)&0xff]^expansion_table0[plaintext_left&0xff]^expansion_table2[(plaintext_left>>16)&0xff]^expansion_table3[plaintext_left>>24]^subkey;
				plaintext_right^=s_permutation_table0[temp_exp&0x3f]^s_permutation_table1[(temp_exp>>6)&0x3f]^s_permutation_table2[(temp_exp>>12)&0x3f]^s_permutation_table3[(temp_exp>>18)&0x3f]^s_permutation_table4[(temp_exp>>24)&0x3f]^s_permutation_table5[(temp_exp>>30)&0x3f]^s_permutation_table6[(temp_exp>>36)&0x3f]^s_permutation_table7[temp_exp>>42];
			}
		}
		if (plaintext_right==ciphertext_left && plaintext_left==ciphertext_right) plaintext_d[0]=key_real;
		key_real+=524288; // Actually THREAD * BLOCK
		key=key_real;
	}
}
void DES_key_schedule(bit64 key, bit64 subkey[16]) {
	key = PC1_table0[key & 0xff] ^ PC1_table1[(key >> 8) & 0xff] ^ PC1_table2[(key >> 16) & 0xff] ^ PC1_table3[(key >> 24) & 0xff] ^ PC1_table4[(key >> 32) & 0xff] ^ PC1_table5[(key >> 40) & 0xff] ^ PC1_table6[(key >> 48) & 0xff] ^ PC1_table7[(key >> 56) & 0xff];
	//Round 1
	key = ((key << 1) & 0xffffffeffffffe) ^ ((key >> 27) & 0x000000010000001);
	subkey[0] = PC2_table0[key & 0xff] ^ PC2_table1[(key >> 8) & 0xff] ^ PC2_table2[(key >> 16) & 0xff] ^ PC2_table3[(key >> 24) & 0xff] ^ PC2_table4[(key >> 32) & 0xff] ^ PC2_table5[(key >> 40) & 0xff] ^ PC2_table6[(key >> 48) & 0xff];
	//Round 2
	key = ((key << 1) & 0xffffffeffffffe) ^ ((key >> 27) & 0x000000010000001);
	subkey[1] = PC2_table0[key & 0xff] ^ PC2_table1[(key >> 8) & 0xff] ^ PC2_table2[(key >> 16) & 0xff] ^ PC2_table3[(key >> 24) & 0xff] ^ PC2_table4[(key >> 32) & 0xff] ^ PC2_table5[(key >> 40) & 0xff] ^ PC2_table6[(key >> 48) & 0xff];
	//Round 3
	key = ((key << 2) & 0xffffffcffffffc) ^ ((key >> 26) & 0x000000030000003);
	subkey[2] = PC2_table0[key & 0xff] ^ PC2_table1[(key >> 8) & 0xff] ^ PC2_table2[(key >> 16) & 0xff] ^ PC2_table3[(key >> 24) & 0xff] ^ PC2_table4[(key >> 32) & 0xff] ^ PC2_table5[(key >> 40) & 0xff] ^ PC2_table6[(key >> 48) & 0xff];
	//Round 4
	key = ((key << 2) & 0xffffffcffffffc) ^ ((key >> 26) & 0x000000030000003);
	subkey[3] = PC2_table0[key & 0xff] ^ PC2_table1[(key >> 8) & 0xff] ^ PC2_table2[(key >> 16) & 0xff] ^ PC2_table3[(key >> 24) & 0xff] ^ PC2_table4[(key >> 32) & 0xff] ^ PC2_table5[(key >> 40) & 0xff] ^ PC2_table6[(key >> 48) & 0xff];
	//Round 5
	key = ((key << 2) & 0xffffffcffffffc) ^ ((key >> 26) & 0x000000030000003);
	subkey[4] = PC2_table0[key & 0xff] ^ PC2_table1[(key >> 8) & 0xff] ^ PC2_table2[(key >> 16) & 0xff] ^ PC2_table3[(key >> 24) & 0xff] ^ PC2_table4[(key >> 32) & 0xff] ^ PC2_table5[(key >> 40) & 0xff] ^ PC2_table6[(key >> 48) & 0xff];
	//Round 6
	key = ((key << 2) & 0xffffffcffffffc) ^ ((key >> 26) & 0x000000030000003);
	subkey[5] = PC2_table0[key & 0xff] ^ PC2_table1[(key >> 8) & 0xff] ^ PC2_table2[(key >> 16) & 0xff] ^ PC2_table3[(key >> 24) & 0xff] ^ PC2_table4[(key >> 32) & 0xff] ^ PC2_table5[(key >> 40) & 0xff] ^ PC2_table6[(key >> 48) & 0xff];
	//Round 7
	key = ((key << 2) & 0xffffffcffffffc) ^ ((key >> 26) & 0x000000030000003);
	subkey[6] = PC2_table0[key & 0xff] ^ PC2_table1[(key >> 8) & 0xff] ^ PC2_table2[(key >> 16) & 0xff] ^ PC2_table3[(key >> 24) & 0xff] ^ PC2_table4[(key >> 32) & 0xff] ^ PC2_table5[(key >> 40) & 0xff] ^ PC2_table6[(key >> 48) & 0xff];
	//Round 8
	key = ((key << 2) & 0xffffffcffffffc) ^ ((key >> 26) & 0x000000030000003);
	subkey[7] = PC2_table0[key & 0xff] ^ PC2_table1[(key >> 8) & 0xff] ^ PC2_table2[(key >> 16) & 0xff] ^ PC2_table3[(key >> 24) & 0xff] ^ PC2_table4[(key >> 32) & 0xff] ^ PC2_table5[(key >> 40) & 0xff] ^ PC2_table6[(key >> 48) & 0xff];
	//Round 9
	key = ((key << 1) & 0xffffffeffffffe) ^ ((key >> 27) & 0x000000010000001);
	subkey[8] = PC2_table0[key & 0xff] ^ PC2_table1[(key >> 8) & 0xff] ^ PC2_table2[(key >> 16) & 0xff] ^ PC2_table3[(key >> 24) & 0xff] ^ PC2_table4[(key >> 32) & 0xff] ^ PC2_table5[(key >> 40) & 0xff] ^ PC2_table6[(key >> 48) & 0xff];
	//Round 10
	key = ((key << 2) & 0xffffffcffffffc) ^ ((key >> 26) & 0x000000030000003);
	subkey[9] = PC2_table0[key & 0xff] ^ PC2_table1[(key >> 8) & 0xff] ^ PC2_table2[(key >> 16) & 0xff] ^ PC2_table3[(key >> 24) & 0xff] ^ PC2_table4[(key >> 32) & 0xff] ^ PC2_table5[(key >> 40) & 0xff] ^ PC2_table6[(key >> 48) & 0xff];
	//Round 11
	key = ((key << 2) & 0xffffffcffffffc) ^ ((key >> 26) & 0x000000030000003);
	subkey[10] = PC2_table0[key & 0xff] ^ PC2_table1[(key >> 8) & 0xff] ^ PC2_table2[(key >> 16) & 0xff] ^ PC2_table3[(key >> 24) & 0xff] ^ PC2_table4[(key >> 32) & 0xff] ^ PC2_table5[(key >> 40) & 0xff] ^ PC2_table6[(key >> 48) & 0xff];
	//Round 12
	key = ((key << 2) & 0xffffffcffffffc) ^ ((key >> 26) & 0x000000030000003);
	subkey[11] = PC2_table0[key & 0xff] ^ PC2_table1[(key >> 8) & 0xff] ^ PC2_table2[(key >> 16) & 0xff] ^ PC2_table3[(key >> 24) & 0xff] ^ PC2_table4[(key >> 32) & 0xff] ^ PC2_table5[(key >> 40) & 0xff] ^ PC2_table6[(key >> 48) & 0xff];
	//Round 13
	key = ((key << 2) & 0xffffffcffffffc) ^ ((key >> 26) & 0x000000030000003);
	subkey[12] = PC2_table0[key & 0xff] ^ PC2_table1[(key >> 8) & 0xff] ^ PC2_table2[(key >> 16) & 0xff] ^ PC2_table3[(key >> 24) & 0xff] ^ PC2_table4[(key >> 32) & 0xff] ^ PC2_table5[(key >> 40) & 0xff] ^ PC2_table6[(key >> 48) & 0xff];
	//Round 14
	key = ((key << 2) & 0xffffffcffffffc) ^ ((key >> 26) & 0x000000030000003);
	subkey[13] = PC2_table0[key & 0xff] ^ PC2_table1[(key >> 8) & 0xff] ^ PC2_table2[(key >> 16) & 0xff] ^ PC2_table3[(key >> 24) & 0xff] ^ PC2_table4[(key >> 32) & 0xff] ^ PC2_table5[(key >> 40) & 0xff] ^ PC2_table6[(key >> 48) & 0xff];
	//Round 15
	key = ((key << 2) & 0xffffffcffffffc) ^ ((key >> 26) & 0x000000030000003);
	subkey[14] = PC2_table0[key & 0xff] ^ PC2_table1[(key >> 8) & 0xff] ^ PC2_table2[(key >> 16) & 0xff] ^ PC2_table3[(key >> 24) & 0xff] ^ PC2_table4[(key >> 32) & 0xff] ^ PC2_table5[(key >> 40) & 0xff] ^ PC2_table6[(key >> 48) & 0xff];
	//Round 16
	key = ((key << 1) & 0xffffffeffffffe) ^ ((key >> 27) & 0x000000010000001);
	subkey[15] = PC2_table0[key & 0xff] ^ PC2_table1[(key >> 8) & 0xff] ^ PC2_table2[(key >> 16) & 0xff] ^ PC2_table3[(key >> 24) & 0xff] ^ PC2_table4[(key >> 32) & 0xff] ^ PC2_table5[(key >> 40) & 0xff] ^ PC2_table6[(key >> 48) & 0xff];
}
void CTR() {
	bit32* expansion_table1_d;
	bit64 key = 0x752978397493cb70, subkey[16];
	bit64* expansion_table0_d, * expansion_table2_d, * expansion_table3_d, * subkey_d;
	bit32* s_permutation_table0_d, * s_permutation_table1_d, * s_permutation_table2_d, * s_permutation_table3_d, * s_permutation_table4_d, * s_permutation_table5_d, * s_permutation_table6_d, * s_permutation_table7_d;
	bit64* ciphertext_d;
	DES_key_schedule(key, subkey);

	cudaMalloc((void**)&expansion_table0_d, 256 * sizeof(bit64));	cudaMemcpy(expansion_table0_d, expansion_table0, 256 * sizeof(bit64), cudaMemcpyHostToDevice);
	cudaMalloc((void**)&expansion_table1_d, 256 * sizeof(bit32));	cudaMemcpy(expansion_table1_d, expansion_table1, 256 * sizeof(bit32), cudaMemcpyHostToDevice);
	cudaMalloc((void**)&expansion_table2_d, 256 * sizeof(bit64));	cudaMemcpy(expansion_table2_d, expansion_table2, 256 * sizeof(bit64), cudaMemcpyHostToDevice);
	cudaMalloc((void**)&expansion_table3_d, 256 * sizeof(bit64));	cudaMemcpy(expansion_table3_d, expansion_table3, 256 * sizeof(bit64), cudaMemcpyHostToDevice);
	cudaMalloc((void**)&s_permutation_table0_d, 64 * sizeof(bit32));	cudaMemcpy(s_permutation_table0_d, s_permutation_table0, 64 * sizeof(bit32), cudaMemcpyHostToDevice);
	cudaMalloc((void**)&s_permutation_table1_d, 64 * sizeof(bit32));	cudaMemcpy(s_permutation_table1_d, s_permutation_table1, 64 * sizeof(bit32), cudaMemcpyHostToDevice);
	cudaMalloc((void**)&s_permutation_table2_d, 64 * sizeof(bit32));	cudaMemcpy(s_permutation_table2_d, s_permutation_table2, 64 * sizeof(bit32), cudaMemcpyHostToDevice);
	cudaMalloc((void**)&s_permutation_table3_d, 64 * sizeof(bit32));	cudaMemcpy(s_permutation_table3_d, s_permutation_table3, 64 * sizeof(bit32), cudaMemcpyHostToDevice);
	cudaMalloc((void**)&s_permutation_table4_d, 64 * sizeof(bit32));	cudaMemcpy(s_permutation_table4_d, s_permutation_table4, 64 * sizeof(bit32), cudaMemcpyHostToDevice);
	cudaMalloc((void**)&s_permutation_table5_d, 64 * sizeof(bit32));	cudaMemcpy(s_permutation_table5_d, s_permutation_table5, 64 * sizeof(bit32), cudaMemcpyHostToDevice);
	cudaMalloc((void**)&s_permutation_table6_d, 64 * sizeof(bit32));	cudaMemcpy(s_permutation_table6_d, s_permutation_table6, 64 * sizeof(bit32), cudaMemcpyHostToDevice);
	cudaMalloc((void**)&s_permutation_table7_d, 64 * sizeof(bit32));	cudaMemcpy(s_permutation_table7_d, s_permutation_table7, 64 * sizeof(bit32), cudaMemcpyHostToDevice);
	cudaMalloc((void**)&subkey_d, 64 * sizeof(bit64));	cudaMemcpy(subkey_d, subkey, 16 * sizeof(bit64), cudaMemcpyHostToDevice);
	cudaMalloc((void**)&ciphertext_d, THREAD * BLOCK * sizeof(bit64));	//cudaMemset(plaintext_d, 0, THREAD * BLOCK * sizeof(bit64));
	StartCounter();
	DES_CTR << <BLOCK, THREAD >> > (subkey_d, ciphertext_d, expansion_table0_d, expansion_table1_d, expansion_table2_d, expansion_table3_d, s_permutation_table0_d, s_permutation_table1_d, s_permutation_table2_d, s_permutation_table3_d, s_permutation_table4_d, s_permutation_table5_d, s_permutation_table6_d, s_permutation_table7_d);
	cudaMemcpy(ciphertext, ciphertext_d, 2 * sizeof(bit64), cudaMemcpyDeviceToHost);

	printf("Time: %lf\n", GetCounter());
	printf("Ciphertext is: %I64x %I64x\n", ciphertext[0], ciphertext[1]);
//	printf("Time: %u seconds\n", clock() / CLOCKS_PER_SEC);
	// Cleanup
	cudaFree(ciphertext_d); cudaFree(subkey_d);
	cudaFree(expansion_table0_d); cudaFree(expansion_table1_d); cudaFree(expansion_table2_d); cudaFree(expansion_table3_d);
	cudaFree(s_permutation_table0_d); cudaFree(s_permutation_table1_d); cudaFree(s_permutation_table2_d); cudaFree(s_permutation_table3_d);
	cudaFree(s_permutation_table4_d); cudaFree(s_permutation_table5_d); cudaFree(s_permutation_table6_d); cudaFree(s_permutation_table7_d);
	printf("%s\n", cudaGetErrorString(cudaGetLastError()));
}
void exhaustive() {
	bit32 *expansion_table1_d,*PC2_table0_d,*PC2_table1_d,*PC2_table2_d,plaintext_left,plaintext_right,ciphertext_left,ciphertext_right;
	bit64 plaintext=0x1122334455667788, ciphertext=0xb5219ee81aa7499d, key=0x752978397493cb70,*plaintext_d, k[2];
	bit64 *expansion_table0_d,*expansion_table2_d,*expansion_table3_d;
	bit64 *PC2_table3_d,*PC2_table4_d,*PC2_table5_d,*PC2_table6_d;
	bit32 *s_permutation_table0_d,*s_permutation_table1_d,*s_permutation_table2_d,*s_permutation_table3_d,*s_permutation_table4_d,*s_permutation_table5_d,*s_permutation_table6_d,*s_permutation_table7_d;
	key=PC1_table0[key&0xff]^PC1_table1[(key>>8)&0xff]^PC1_table2[(key>>16)&0xff]^PC1_table3[(key>>24)&0xff]^PC1_table4[(key>>32)&0xff]^PC1_table5[(key>>40)&0xff]^PC1_table6[(key>>48)&0xff]^PC1_table7[(key>>56)&0xff];
	plaintext=IP_table0[plaintext&0xff]^IP_table1[(plaintext>>8)&0xff]^IP_table2[(plaintext>>16)&0xff]^IP_table3[(plaintext>>24)&0xff]^IP_table4[(plaintext>>32)&0xff]^IP_table5[(plaintext>>40)&0xff]^IP_table6[(plaintext>>48)&0xff]^IP_table7[(plaintext>>56)&0xff];
	ciphertext=FP2_table0[ciphertext&0xff]^FP2_table1[(ciphertext>>8)&0xff]^FP2_table2[(ciphertext>>16)&0xff]^FP2_table3[(ciphertext>>24)&0xff]^FP2_table4[(ciphertext>>32)&0xff]^FP2_table5[(ciphertext>>40)&0xff]^FP2_table6[(ciphertext>>48)&0xff]^FP2_table7[(ciphertext>>56)&0xff];
	plaintext_left=plaintext>>32; plaintext_right=plaintext&0xffffffff;
	ciphertext_left=ciphertext>>32; ciphertext_right=ciphertext&0xffffffff;
	cudaMalloc((void **)&expansion_table0_d, 256*sizeof(bit64));	cudaMemcpy(expansion_table0_d,expansion_table0,256*sizeof(bit64),cudaMemcpyHostToDevice);
	cudaMalloc((void **)&expansion_table1_d, 256*sizeof(bit32));	cudaMemcpy(expansion_table1_d,expansion_table1,256*sizeof(bit32),cudaMemcpyHostToDevice);
	cudaMalloc((void **)&expansion_table2_d, 256*sizeof(bit64));	cudaMemcpy(expansion_table2_d,expansion_table2,256*sizeof(bit64),cudaMemcpyHostToDevice);
	cudaMalloc((void **)&expansion_table3_d, 256*sizeof(bit64));	cudaMemcpy(expansion_table3_d,expansion_table3,256*sizeof(bit64),cudaMemcpyHostToDevice);
	cudaMalloc((void **)&PC2_table0_d, 256*sizeof(bit32));			cudaMemcpy(PC2_table0_d,PC2_table0,256*sizeof(bit32),cudaMemcpyHostToDevice);
	cudaMalloc((void **)&PC2_table1_d, 256*sizeof(bit32));			cudaMemcpy(PC2_table1_d,PC2_table1,256*sizeof(bit32),cudaMemcpyHostToDevice);
	cudaMalloc((void **)&PC2_table2_d, 256*sizeof(bit32));			cudaMemcpy(PC2_table2_d,PC2_table2,256*sizeof(bit32),cudaMemcpyHostToDevice);
	cudaMalloc((void **)&PC2_table3_d, 256*sizeof(bit64));			cudaMemcpy(PC2_table3_d,PC2_table3,256*sizeof(bit64),cudaMemcpyHostToDevice);
	cudaMalloc((void **)&PC2_table4_d, 256*sizeof(bit64));			cudaMemcpy(PC2_table4_d,PC2_table4,256*sizeof(bit64),cudaMemcpyHostToDevice);
	cudaMalloc((void **)&PC2_table5_d, 256*sizeof(bit64));			cudaMemcpy(PC2_table5_d,PC2_table5,256*sizeof(bit64),cudaMemcpyHostToDevice);
	cudaMalloc((void **)&PC2_table6_d, 256*sizeof(bit64));			cudaMemcpy(PC2_table6_d,PC2_table6,256*sizeof(bit64),cudaMemcpyHostToDevice);
	cudaMalloc((void **)&s_permutation_table0_d, 64*sizeof(bit32));	cudaMemcpy(s_permutation_table0_d,s_permutation_table0,64*sizeof(bit32),cudaMemcpyHostToDevice);
	cudaMalloc((void **)&s_permutation_table1_d, 64*sizeof(bit32));	cudaMemcpy(s_permutation_table1_d,s_permutation_table1,64*sizeof(bit32),cudaMemcpyHostToDevice);
	cudaMalloc((void **)&s_permutation_table2_d, 64*sizeof(bit32));	cudaMemcpy(s_permutation_table2_d,s_permutation_table2,64*sizeof(bit32),cudaMemcpyHostToDevice);
	cudaMalloc((void **)&s_permutation_table3_d, 64*sizeof(bit32));	cudaMemcpy(s_permutation_table3_d,s_permutation_table3,64*sizeof(bit32),cudaMemcpyHostToDevice);
	cudaMalloc((void **)&s_permutation_table4_d, 64*sizeof(bit32));	cudaMemcpy(s_permutation_table4_d,s_permutation_table4,64*sizeof(bit32),cudaMemcpyHostToDevice);
	cudaMalloc((void **)&s_permutation_table5_d, 64*sizeof(bit32));	cudaMemcpy(s_permutation_table5_d,s_permutation_table5,64*sizeof(bit32),cudaMemcpyHostToDevice);
	cudaMalloc((void **)&s_permutation_table6_d, 64*sizeof(bit32));	cudaMemcpy(s_permutation_table6_d,s_permutation_table6,64*sizeof(bit32),cudaMemcpyHostToDevice);
	cudaMalloc((void **)&s_permutation_table7_d, 64*sizeof(bit32));	cudaMemcpy(s_permutation_table7_d,s_permutation_table7,64*sizeof(bit32),cudaMemcpyHostToDevice);
	cudaMalloc((void **)&plaintext_d, 2*sizeof(bit64));	cudaMemset(plaintext_d,0,2*sizeof(bit64));
	StartCounter();
	DES_exhaustive<<<BLOCK, THREAD >>>(plaintext_d,key,plaintext_left,plaintext_right,ciphertext_left,ciphertext_right,expansion_table0_d,expansion_table1_d,expansion_table2_d,expansion_table3_d,PC2_table0_d,PC2_table1_d,PC2_table2_d,PC2_table3_d,PC2_table4_d,PC2_table5_d,PC2_table6_d,s_permutation_table0_d,s_permutation_table1_d,s_permutation_table2_d,s_permutation_table3_d,s_permutation_table4_d,s_permutation_table5_d,s_permutation_table6_d,s_permutation_table7_d);
	cudaMemcpy(k,plaintext_d,2*sizeof(bit64),cudaMemcpyDeviceToHost);
	if (k[0] || k[1]) printf("Correct key is: %I64x %I64x\n",k[1], k[0]);
	printf("Time: %lf\n", GetCounter());
	printf("Time: %u seconds\n", clock() / CLOCKS_PER_SEC);
	// Cleanup
	cudaFree(plaintext_d);cudaFree(expansion_table0_d);cudaFree(expansion_table1_d);cudaFree(expansion_table2_d);cudaFree(expansion_table3_d);
	cudaFree(PC2_table0_d);cudaFree(PC2_table1_d);cudaFree(PC2_table2_d);cudaFree(PC2_table3_d);cudaFree(PC2_table4_d);cudaFree(PC2_table5_d);cudaFree(PC2_table6_d);
	cudaFree(s_permutation_table0_d);cudaFree(s_permutation_table1_d);cudaFree(s_permutation_table2_d);cudaFree(s_permutation_table3_d);
	cudaFree(s_permutation_table4_d);cudaFree(s_permutation_table5_d);cudaFree(s_permutation_table6_d);cudaFree(s_permutation_table7_d);
	printf("%s\n",cudaGetErrorString(cudaGetLastError()));
}

int main() {
	int choice = 0;
	printf(
		"(1) Exhaustive search\n"
		"(2) Counter mode\n"
		"Enter choice: "
	);
	scanf_s("%d", &choice);
	if (choice == 1) exhaustive();
	if (choice == 2) CTR();
	return 1;
}
