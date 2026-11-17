Characteristic:				Baseline A		Baseline B
Modulation:				8-PSK Gray		8-PSK Gray
Turbo Encoder:				No			Yes
Channel:				AWGN			AWGN
Soft Demapping:				Not necessary		Yes
Turbo Decoder:				No			Yes
Result:					BER			BER and PER
External Reference:			berawgn			Not present yet
Quantitative validation with theory:	Yes			Not yet

-------------------------------------------------------------------------------------------------

Characteristic:		Baseline B				Scientifical Article
Modulation:		8-PSK Gray				8-PSK con Gray and other mappings
Channel:		AWGN					AWGN
Turbo Coding:		Nominal rate 1/3			Rate 1/3
Info length:		1331 bits				1530 bits
RSC Coders:		Default MATLAB				16 states, polynoms (23,31)
Turbo Iterations:	Default MATLAB				8
Demapping:		Soft, LLR approximate			Iterative demapping BICM-ID
Interleaver:		Casual per packet			Pseudocasual
Objective:		BER/PER of the coded chain 		BER comparison between different mappings