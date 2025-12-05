`timescale 1ns/1ps
module bch(
	input clk,
	input rstn,
	input mode,
	input [1:0] code,
	input set,
	input [63:0] idata,
	output ready,
	output finish,
	output [9:0] odata
);

	localparam S_IDLE = 0;
	localparam S_LOAD = 1;
	localparam S_SYN  = 2;
	localparam S_BER  = 3;
	localparam S_CHI_HARD  = 4;
	localparam S_CHI_SOFT1 = 5;
	localparam S_CHI_SOFT2 = 6;
	localparam S_CORR = 7;
	localparam S_OUT  = 8;

	localparam ALPHA_1 = 11'b00000000010;
	localparam ALPHA_2 = 11'b00000000100;
	localparam ALPHA_3 = 11'b00000001000;
	localparam ALPHA_4 = 11'b00000010000;
	localparam ALPHA_5 = 11'b00000100000;
	localparam ALPHA_6 = 11'b00001000000;
	localparam ALPHA_7 = 11'b00010000000;
	localparam ALPHA_8 = 11'b00100000000;

    localparam ALPHA_NEG1_6  = 11'b00000100001;
    localparam ALPHA_NEG2_6  = 11'b00000110001;
    localparam ALPHA_NEG3_6  = 11'b00000111001;
    localparam ALPHA_NEG4_6  = 11'b00000111101;
    localparam ALPHA_NEG5_6  = 11'b00000111111;
    localparam ALPHA_NEG6_6  = 11'b00000111110;
    localparam ALPHA_NEG7_6  = 11'b00000011111;
    localparam ALPHA_NEG8_6  = 11'b00000101110;
	localparam ALPHA_NEG16_6 = 11'b00000100111;
	localparam ALPHA_NEG24_6 = 11'b00000110110;
	localparam ALPHA_NEG32_6 = 11'b00000100101;

	localparam ALPHA_NEG1_8  = 11'b00010001110;
	localparam ALPHA_NEG2_8  = 11'b00001000111;
    localparam ALPHA_NEG3_8  = 11'b00010101101;
    localparam ALPHA_NEG4_8  = 11'b00011011000;
    localparam ALPHA_NEG5_8  = 11'b00001101100;
    localparam ALPHA_NEG6_8  = 11'b00000110110;
    localparam ALPHA_NEG7_8  = 11'b00000011011;
	localparam ALPHA_NEG8_8  = 11'b00010000011;
	localparam ALPHA_NEG16_8 = 11'b00000010110;
	localparam ALPHA_NEG24_8 = 11'b00011110101;
	localparam ALPHA_NEG32_8 = 11'b00000001001;

	localparam ALPHA_NEG1_10  = 11'b01000000100;
	localparam ALPHA_NEG2_10  = 11'b00100000010;
    localparam ALPHA_NEG3_10  = 11'b00010000001;
    localparam ALPHA_NEG4_10  = 11'b01001000100;
    localparam ALPHA_NEG5_10  = 11'b00100100010;
    localparam ALPHA_NEG6_10  = 11'b00010010001;
    localparam ALPHA_NEG7_10  = 11'b01001001100;
	localparam ALPHA_NEG8_10  = 11'b00100100110;
	localparam ALPHA_NEG16_10 = 11'b01001011101;
	localparam ALPHA_NEG24_10 = 11'b01111010001;
	localparam ALPHA_NEG32_10 = 11'b00001100111;

	reg ready_r, ready_w;

	reg [3:0] state_r, state_w;
	reg [9:0] cnt_r, cnt_w;

	// input data
	reg [7:0] data_r [0:1023], data_w [0:1023];
	reg mode_r, mode_w;
	reg [1:0] code_r, code_w;

	// syndrome calculation
	reg [10:0] S_r [0:7], S_w [0:7];
	wire S1_S2_0 = (S_r[0] == 0 && S_r[1] == 0);
	wire S3_S4_0 = (S_r[2] == 0 && S_r[3] == 0);
	wire S5_S6_0 = (S_r[4] == 0 && S_r[5] == 0);
	wire S7_S8_0 = (S_r[6] == 0 && S_r[7] == 0);
	wire S1_S4_0 = S1_S2_0 && S3_S4_0;
	wire S5_S8_0 = S5_S6_0 && S7_S8_0;
	wire S1_S8_0 = S1_S4_0 && S5_S8_0;
	reg [9:0] S_temp1_w [0:7], S_temp2_w [0:7], S_temp3_w [0:7], S_temp4_w [0:7], S_temp5_w [0:7], S_temp6_w [0:7], S_temp7_w [0:7];
	
	// soft decision
	wire [7:0] abs_idata [0:7];
	assign abs_idata[0] = idata[7] ? (~idata[7:0] + 1) : idata[7:0];
	assign abs_idata[1] = idata[15] ? (~idata[15:8] + 1) : idata[15:8];
	assign abs_idata[2] = idata[23] ? (~idata[23:16] + 1) : idata[23:16];
	assign abs_idata[3] = idata[31] ? (~idata[31:24] + 1) : idata[31:24];
	assign abs_idata[4] = idata[39] ? (~idata[39:32] + 1) : idata[39:32];
	assign abs_idata[5] = idata[47] ? (~idata[47:40] + 1) : idata[47:40];
	assign abs_idata[6] = idata[55] ? (~idata[55:48] + 1) : idata[55:48];
	assign abs_idata[7] = idata[63] ? (~idata[63:56] + 1) : idata[63:56];
	wire [7:0] abs_data7 = (((code_r == 3) && (cnt_r == 1023)) || ((code_r == 2) && (cnt_r == 255)) || ((code_r == 1) && (cnt_r == 63))) ? 8'b11111111 : abs_idata[7];
	reg [7:0] min1_r, min2_r, min1_w, min2_w, min1_out, min2_out;
	reg [9:0] index1_r, index2_r, index1_w, index2_w, index1_out, index2_out;
	wire [63:0] abs_data = {abs_data7, abs_idata[6], abs_idata[5], abs_idata[4], abs_idata[3], abs_idata[2], abs_idata[1], abs_idata[0]};

	find_min Min(
	.abs_data(abs_data), .cnt(cnt_r), .min1_in(min1_r), .min2_in(min2_r), .index1_in(index1_r), .index2_in(index2_r),
	.min1_out(min1_out), .min2_out(min2_out), .index1_out(index1_out), .index2_out(index2_out)
);

	//  ber algorithm
	reg [10:0] delta_r [0:4], delta_w [0:4], delta_rho_r [0:4], delta_rho_w [0:4], temp1_w [0:4], temp2_w [0:4], temp3_w[0:4];
	reg [10:0] d_r, d_w, d_rho_r, d_rho_w; 
	reg [3:0]  l_r, l_w, l_rho_r, l_rho_w;
	reg [3:0]  rho_r, rho_w;
	reg [3:0]  cnt_temp_w;

	// chien search
	// reg [10:0] power_r [0:7], power_w [0:7];
	// reg [10:0] power2_r [0:7], power2_w [0:7];
	// reg [10:0] power3_r [0:7], power3_w [0:7];
	// reg [10:0] power4_r [0:7], power4_w [0:7];
	reg [10:0] stage1_buff_r [0:7][0:4], stage1_buff_w [0:7][0:4]; // first index: parallel number
	reg [2:0]  root_cnt_r, root_cnt_w;
	reg [10:0] temp_root_w [0:7];
	reg [2:0]  temp_root_cnt_w [0:7];
	reg [9:0]  root_r [0:3], root_w [0:3];
	wire [10:0] delta_poly_w [0:7][0:4];    // first index: parallel number, second index: degree

	genvar k;
	generate
		for (k = 0; k < 8; k = k + 1) begin
			assign delta_poly_w[k][0] = delta_r[0];
		end
		for (k = 1; k < 5; k = k + 1) begin
			assign delta_poly_w[0][k] = delta_r[k];
		end
		assign delta_poly_w[1][1] = shift_poly_1(delta_r[1]);
		assign delta_poly_w[1][2] = shift_poly_2(delta_r[2]);
		assign delta_poly_w[1][3] = shift_poly_3(delta_r[3]);
		assign delta_poly_w[1][4] = shift_poly_4(delta_r[4]);
		assign delta_poly_w[2][1] = shift_poly_2(delta_r[1]);
		assign delta_poly_w[2][2] = shift_poly_4(delta_r[2]);
		assign delta_poly_w[2][3] = shift_poly_6(delta_r[3]);
		assign delta_poly_w[2][4] = shift_poly_8(delta_r[4]);
		assign delta_poly_w[3][1] = shift_poly_3(delta_r[1]);
		assign delta_poly_w[3][2] = shift_poly_6(delta_r[2]);
		assign delta_poly_w[3][3] = shift_poly_4(shift_poly_5(delta_r[3]));
		assign delta_poly_w[3][4] = shift_poly_6(shift_poly_6(delta_r[4]));
		assign delta_poly_w[4][1] = shift_poly_4(delta_r[1]);
		assign delta_poly_w[4][2] = shift_poly_8(delta_r[2]);
		assign delta_poly_w[4][3] = shift_poly_6(shift_poly_6(delta_r[3]));
		assign delta_poly_w[4][4] = shift_poly_8(shift_poly_8(delta_r[4]));
		assign delta_poly_w[5][1] = shift_poly_5(delta_r[1]);
		assign delta_poly_w[5][2] = shift_poly_5(shift_poly_5(delta_r[2]));
		assign delta_poly_w[5][3] = shift_poly_7(shift_poly_8(delta_r[3]));
		assign delta_poly_w[5][4] = shift_poly_6(shift_poly_7(shift_poly_7(delta_r[4])));
		assign delta_poly_w[6][1] = shift_poly_6(delta_r[1]);
		assign delta_poly_w[6][2] = shift_poly_6(shift_poly_6(delta_r[2]));
		assign delta_poly_w[6][3] = shift_poly_6(shift_poly_6(shift_poly_6(delta_r[3])));
		assign delta_poly_w[6][4] = shift_poly_8(shift_poly_8(shift_poly_8(delta_r[4])));
		assign delta_poly_w[7][1] = shift_poly_7(delta_r[1]);
		assign delta_poly_w[7][2] = shift_poly_7(shift_poly_7(delta_r[2]));
		assign delta_poly_w[7][3] = shift_poly_7(shift_poly_7(shift_poly_7(delta_r[3])));
		assign delta_poly_w[7][4] = shift_poly_7(shift_poly_7(shift_poly_7(shift_poly_7(delta_r[4]))));
	endgenerate

	// correlation
	reg [6:0] corr_r[0:3], corr_w[0:3]; // 0 for first stage and 1 - 3 for second stage

	// output 
	reg [9:0] odata_r, odata_w;
	reg finish_r, finish_w;
	assign odata = odata_r;
	assign finish = finish_r;

	integer i, j;

	// clock gating
	wire load_en     = 1;
	wire syn_low_en  = 1;
	wire syn_high_en = 1;
	wire ber_en      = 1;

	assign ready = ready_r;

	always @(*) begin
		state_w = state_r;
		cnt_w = cnt_r;
		ready_w = ready_r;
		code_w = code_r;
		mode_w = mode_r;
		d_w = d_r; 
		d_rho_w = d_rho_r;
		l_w = l_r; 
		l_rho_w = l_rho_r;
		rho_w = rho_r;
		finish_w = finish_r;;
		odata_w = odata_r;
		min1_w = min1_r;
		min2_w = min2_r;
		index1_w = index1_r;
		index2_w = index2_r;
		cnt_temp_w = 0;
		for (i = 0; i < 8; i = i + 1) S_w[i] = S_r[i];
		for (i = 0; i < 8; i = i + 1) S_temp1_w[i] = 0;
		for (i = 0; i < 8; i = i + 1) S_temp2_w[i] = 0;
		for (i = 0; i < 8; i = i + 1) S_temp3_w[i] = 0;
		for (i = 0; i < 8; i = i + 1) S_temp4_w[i] = 0;
		for (i = 0; i < 8; i = i + 1) S_temp5_w[i] = 0;
		for (i = 0; i < 8; i = i + 1) S_temp6_w[i] = 0;
		for (i = 0; i < 8; i = i + 1) S_temp7_w[i] = 0;
		for (i = 0; i < 5; i = i + 1) begin
			delta_w[i] = delta_r[i];
			delta_rho_w[i] = delta_rho_r[i];
			temp1_w[i] = 11'b0;
			temp2_w[i] = 11'b0;
			temp3_w[i] = 11'b0;
		end
		for (i = 0; i < 1024; i = i + 1) begin
			data_w[i] = data_r[i];
		end 
		for (i = 0; i < 8; i = i + 1) begin
			S_w[i] = S_r[i];
			temp_root_w[i] = 1;
			temp_root_cnt_w[i] = 0;
			// power_w[i] = power_r[i];
			// power2_w[i] = power2_r[i];
			// power3_w[i] = power3_r[i];
			// power4_w[i] = power4_r[i];
			for (j = 0; j < 5; j = j + 1) begin
				stage1_buff_w[i][j] = stage1_buff_r[i][j];
			end
		end 
		for (i = 0; i < 4; i = i + 1) begin
			root_w[i] = root_r[i];
			corr_w[i] = corr_r[i];
		end
		root_cnt_w = root_cnt_r;
		case (state_r)
			S_IDLE: begin
				if (set && !ready_r) begin
					case (code)
						1: cnt_w = 63;
						2: cnt_w = 255;
						3: cnt_w = 1023;
						default: cnt_w = 1023;
					endcase
					code_w = code;
					mode_w = mode;
					ready_w = 1'b1;
					state_w = S_LOAD;
				end
				// else if (ready_r) state_w = S_LOAD;
				else state_w = S_IDLE;

				min1_w = 8'b11111111;
				min2_w = 8'b11111111;

			end  
			S_LOAD: begin
				case (code_r)
					1: begin
						if (cnt_r == 63) begin
							S_temp2_w[0] = {9'b0, idata[55]};
							S_temp2_w[1] = {9'b0, idata[55]};
							S_temp2_w[2] = {9'b0, idata[55]};
							S_temp2_w[3] = {9'b0, idata[55]};
							S_temp3_w[0] = shift_poly_1(S_temp2_w[0]) ^ {9'b0, idata[47]};
							S_temp3_w[1] = shift_poly_2(S_temp2_w[1]) ^ {9'b0, idata[47]};
							S_temp3_w[2] = shift_poly_3(S_temp2_w[2]) ^ {9'b0, idata[47]};
							S_temp3_w[3] = shift_poly_4(S_temp2_w[3]) ^ {9'b0, idata[47]};
							S_temp4_w[0] = shift_poly_1(S_temp3_w[0]) ^ {9'b0, idata[39]};
							S_temp4_w[1] = shift_poly_2(S_temp3_w[1]) ^ {9'b0, idata[39]};
							S_temp4_w[2] = shift_poly_3(S_temp3_w[2]) ^ {9'b0, idata[39]};
							S_temp4_w[3] = shift_poly_4(S_temp3_w[3]) ^ {9'b0, idata[39]};
							S_temp5_w[0] = shift_poly_1(S_temp4_w[0]) ^ {9'b0, idata[31]};
							S_temp5_w[1] = shift_poly_2(S_temp4_w[1]) ^ {9'b0, idata[31]};
							S_temp5_w[2] = shift_poly_3(S_temp4_w[2]) ^ {9'b0, idata[31]};
							S_temp5_w[3] = shift_poly_4(S_temp4_w[3]) ^ {9'b0, idata[31]};
							S_temp6_w[0] = shift_poly_1(S_temp5_w[0]) ^ {9'b0, idata[23]};
							S_temp6_w[1] = shift_poly_2(S_temp5_w[1]) ^ {9'b0, idata[23]};
							S_temp6_w[2] = shift_poly_3(S_temp5_w[2]) ^ {9'b0, idata[23]};
							S_temp6_w[3] = shift_poly_4(S_temp5_w[3]) ^ {9'b0, idata[23]};
							S_temp7_w[0] = shift_poly_1(S_temp6_w[0]) ^ {9'b0, idata[15]};
							S_temp7_w[1] = shift_poly_2(S_temp6_w[1]) ^ {9'b0, idata[15]};
							S_temp7_w[2] = shift_poly_3(S_temp6_w[2]) ^ {9'b0, idata[15]};
							S_temp7_w[3] = shift_poly_4(S_temp6_w[3]) ^ {9'b0, idata[15]};
							S_w[0] = shift_poly_1(S_temp7_w[0]) ^ {9'b0, idata[7]};
							S_w[1] = shift_poly_2(S_temp7_w[1]) ^ {9'b0, idata[7]};
							S_w[2] = shift_poly_3(S_temp7_w[2]) ^ {9'b0, idata[7]};
							S_w[3] = shift_poly_4(S_temp7_w[3]) ^ {9'b0, idata[7]};
						end
						else if (cnt_r >= 7) begin
							S_temp1_w[0] = shift_poly_1(S_r[0]) ^ {9'b0, idata[63]};
							S_temp1_w[1] = shift_poly_2(S_r[1]) ^ {9'b0, idata[63]};
							S_temp1_w[2] = shift_poly_3(S_r[2]) ^ {9'b0, idata[63]};
							S_temp1_w[3] = shift_poly_4(S_r[3]) ^ {9'b0, idata[63]};
							S_temp2_w[0] = shift_poly_1(S_temp1_w[0]) ^ {9'b0, idata[55]};
							S_temp2_w[1] = shift_poly_2(S_temp1_w[1]) ^ {9'b0, idata[55]};
							S_temp2_w[2] = shift_poly_3(S_temp1_w[2]) ^ {9'b0, idata[55]};
							S_temp2_w[3] = shift_poly_4(S_temp1_w[3]) ^ {9'b0, idata[55]};
							S_temp3_w[0] = shift_poly_1(S_temp2_w[0]) ^ {9'b0, idata[47]};
							S_temp3_w[1] = shift_poly_2(S_temp2_w[1]) ^ {9'b0, idata[47]};
							S_temp3_w[2] = shift_poly_3(S_temp2_w[2]) ^ {9'b0, idata[47]};
							S_temp3_w[3] = shift_poly_4(S_temp2_w[3]) ^ {9'b0, idata[47]};
							S_temp4_w[0] = shift_poly_1(S_temp3_w[0]) ^ {9'b0, idata[39]};
							S_temp4_w[1] = shift_poly_2(S_temp3_w[1]) ^ {9'b0, idata[39]};
							S_temp4_w[2] = shift_poly_3(S_temp3_w[2]) ^ {9'b0, idata[39]};
							S_temp4_w[3] = shift_poly_4(S_temp3_w[3]) ^ {9'b0, idata[39]};
							S_temp5_w[0] = shift_poly_1(S_temp4_w[0]) ^ {9'b0, idata[31]};
							S_temp5_w[1] = shift_poly_2(S_temp4_w[1]) ^ {9'b0, idata[31]};
							S_temp5_w[2] = shift_poly_3(S_temp4_w[2]) ^ {9'b0, idata[31]};
							S_temp5_w[3] = shift_poly_4(S_temp4_w[3]) ^ {9'b0, idata[31]};
							S_temp6_w[0] = shift_poly_1(S_temp5_w[0]) ^ {9'b0, idata[23]};
							S_temp6_w[1] = shift_poly_2(S_temp5_w[1]) ^ {9'b0, idata[23]};
							S_temp6_w[2] = shift_poly_3(S_temp5_w[2]) ^ {9'b0, idata[23]};
							S_temp6_w[3] = shift_poly_4(S_temp5_w[3]) ^ {9'b0, idata[23]};
							S_temp7_w[0] = shift_poly_1(S_temp6_w[0]) ^ {9'b0, idata[15]};
							S_temp7_w[1] = shift_poly_2(S_temp6_w[1]) ^ {9'b0, idata[15]};
							S_temp7_w[2] = shift_poly_3(S_temp6_w[2]) ^ {9'b0, idata[15]};
							S_temp7_w[3] = shift_poly_4(S_temp6_w[3]) ^ {9'b0, idata[15]};
							S_w[0] = shift_poly_1(S_temp7_w[0]) ^ {9'b0, idata[7]};
							S_w[1] = shift_poly_2(S_temp7_w[1]) ^ {9'b0, idata[7]};
							S_w[2] = shift_poly_3(S_temp7_w[2]) ^ {9'b0, idata[7]};
							S_w[3] = shift_poly_4(S_temp7_w[3]) ^ {9'b0, idata[7]};
						end
					end
					2: begin
						if (cnt_r == 255) begin
							S_temp2_w[0] = {9'b0, idata[55]};
							S_temp2_w[1] = {9'b0, idata[55]};
							S_temp2_w[2] = {9'b0, idata[55]};
							S_temp2_w[3] = {9'b0, idata[55]};
							S_temp3_w[0] = shift_poly_1(S_temp2_w[0]) ^ {9'b0, idata[47]};
							S_temp3_w[1] = shift_poly_2(S_temp2_w[1]) ^ {9'b0, idata[47]};
							S_temp3_w[2] = shift_poly_3(S_temp2_w[2]) ^ {9'b0, idata[47]};
							S_temp3_w[3] = shift_poly_4(S_temp2_w[3]) ^ {9'b0, idata[47]};
							S_temp4_w[0] = shift_poly_1(S_temp3_w[0]) ^ {9'b0, idata[39]};
							S_temp4_w[1] = shift_poly_2(S_temp3_w[1]) ^ {9'b0, idata[39]};
							S_temp4_w[2] = shift_poly_3(S_temp3_w[2]) ^ {9'b0, idata[39]};
							S_temp4_w[3] = shift_poly_4(S_temp3_w[3]) ^ {9'b0, idata[39]};
							S_temp5_w[0] = shift_poly_1(S_temp4_w[0]) ^ {9'b0, idata[31]};
							S_temp5_w[1] = shift_poly_2(S_temp4_w[1]) ^ {9'b0, idata[31]};
							S_temp5_w[2] = shift_poly_3(S_temp4_w[2]) ^ {9'b0, idata[31]};
							S_temp5_w[3] = shift_poly_4(S_temp4_w[3]) ^ {9'b0, idata[31]};
							S_temp6_w[0] = shift_poly_1(S_temp5_w[0]) ^ {9'b0, idata[23]};
							S_temp6_w[1] = shift_poly_2(S_temp5_w[1]) ^ {9'b0, idata[23]};
							S_temp6_w[2] = shift_poly_3(S_temp5_w[2]) ^ {9'b0, idata[23]};
							S_temp6_w[3] = shift_poly_4(S_temp5_w[3]) ^ {9'b0, idata[23]};
							S_temp7_w[0] = shift_poly_1(S_temp6_w[0]) ^ {9'b0, idata[15]};
							S_temp7_w[1] = shift_poly_2(S_temp6_w[1]) ^ {9'b0, idata[15]};
							S_temp7_w[2] = shift_poly_3(S_temp6_w[2]) ^ {9'b0, idata[15]};
							S_temp7_w[3] = shift_poly_4(S_temp6_w[3]) ^ {9'b0, idata[15]};
							S_w[0] = shift_poly_1(S_temp7_w[0]) ^ {9'b0, idata[7]};
							S_w[1] = shift_poly_2(S_temp7_w[1]) ^ {9'b0, idata[7]};
							S_w[2] = shift_poly_3(S_temp7_w[2]) ^ {9'b0, idata[7]};
							S_w[3] = shift_poly_4(S_temp7_w[3]) ^ {9'b0, idata[7]};
						end
						else if (cnt_r >= 7) begin
							S_temp1_w[0] = shift_poly_1(S_r[0]) ^ {9'b0, idata[63]};
							S_temp1_w[1] = shift_poly_2(S_r[1]) ^ {9'b0, idata[63]};
							S_temp1_w[2] = shift_poly_3(S_r[2]) ^ {9'b0, idata[63]};
							S_temp1_w[3] = shift_poly_4(S_r[3]) ^ {9'b0, idata[63]};
							S_temp2_w[0] = shift_poly_1(S_temp1_w[0]) ^ {9'b0, idata[55]};
							S_temp2_w[1] = shift_poly_2(S_temp1_w[1]) ^ {9'b0, idata[55]};
							S_temp2_w[2] = shift_poly_3(S_temp1_w[2]) ^ {9'b0, idata[55]};
							S_temp2_w[3] = shift_poly_4(S_temp1_w[3]) ^ {9'b0, idata[55]};
							S_temp3_w[0] = shift_poly_1(S_temp2_w[0]) ^ {9'b0, idata[47]};
							S_temp3_w[1] = shift_poly_2(S_temp2_w[1]) ^ {9'b0, idata[47]};
							S_temp3_w[2] = shift_poly_3(S_temp2_w[2]) ^ {9'b0, idata[47]};
							S_temp3_w[3] = shift_poly_4(S_temp2_w[3]) ^ {9'b0, idata[47]};
							S_temp4_w[0] = shift_poly_1(S_temp3_w[0]) ^ {9'b0, idata[39]};
							S_temp4_w[1] = shift_poly_2(S_temp3_w[1]) ^ {9'b0, idata[39]};
							S_temp4_w[2] = shift_poly_3(S_temp3_w[2]) ^ {9'b0, idata[39]};
							S_temp4_w[3] = shift_poly_4(S_temp3_w[3]) ^ {9'b0, idata[39]};
							S_temp5_w[0] = shift_poly_1(S_temp4_w[0]) ^ {9'b0, idata[31]};
							S_temp5_w[1] = shift_poly_2(S_temp4_w[1]) ^ {9'b0, idata[31]};
							S_temp5_w[2] = shift_poly_3(S_temp4_w[2]) ^ {9'b0, idata[31]};
							S_temp5_w[3] = shift_poly_4(S_temp4_w[3]) ^ {9'b0, idata[31]};
							S_temp6_w[0] = shift_poly_1(S_temp5_w[0]) ^ {9'b0, idata[23]};
							S_temp6_w[1] = shift_poly_2(S_temp5_w[1]) ^ {9'b0, idata[23]};
							S_temp6_w[2] = shift_poly_3(S_temp5_w[2]) ^ {9'b0, idata[23]};
							S_temp6_w[3] = shift_poly_4(S_temp5_w[3]) ^ {9'b0, idata[23]};
							S_temp7_w[0] = shift_poly_1(S_temp6_w[0]) ^ {9'b0, idata[15]};
							S_temp7_w[1] = shift_poly_2(S_temp6_w[1]) ^ {9'b0, idata[15]};
							S_temp7_w[2] = shift_poly_3(S_temp6_w[2]) ^ {9'b0, idata[15]};
							S_temp7_w[3] = shift_poly_4(S_temp6_w[3]) ^ {9'b0, idata[15]};
							S_w[0] = shift_poly_1(S_temp7_w[0]) ^ {9'b0, idata[7]};
							S_w[1] = shift_poly_2(S_temp7_w[1]) ^ {9'b0, idata[7]};
							S_w[2] = shift_poly_3(S_temp7_w[2]) ^ {9'b0, idata[7]};
							S_w[3] = shift_poly_4(S_temp7_w[3]) ^ {9'b0, idata[7]};
						end
					end
					3: begin
						if (cnt_r == 1023) begin
							S_temp2_w[0] = {9'b0, idata[55]};
							S_temp2_w[1] = {9'b0, idata[55]};
							S_temp2_w[2] = {9'b0, idata[55]};
							S_temp2_w[3] = {9'b0, idata[55]};
							S_temp2_w[4] = {9'b0, idata[55]};
							S_temp2_w[5] = {9'b0, idata[55]};
							S_temp2_w[6] = {9'b0, idata[55]};
							S_temp2_w[7] = {9'b0, idata[55]};
							S_temp3_w[0] = shift_poly_1(S_temp2_w[0]) ^ {9'b0, idata[47]};
							S_temp3_w[1] = shift_poly_2(S_temp2_w[1]) ^ {9'b0, idata[47]};
							S_temp3_w[2] = shift_poly_3(S_temp2_w[2]) ^ {9'b0, idata[47]};
							S_temp3_w[3] = shift_poly_4(S_temp2_w[3]) ^ {9'b0, idata[47]};
							S_temp3_w[4] = shift_poly_5(S_temp2_w[4]) ^ {9'b0, idata[47]};
							S_temp3_w[5] = shift_poly_6(S_temp2_w[5]) ^ {9'b0, idata[47]};
							S_temp3_w[6] = shift_poly_7(S_temp2_w[6]) ^ {9'b0, idata[47]};
							S_temp3_w[7] = shift_poly_8(S_temp2_w[7]) ^ {9'b0, idata[47]};
							S_temp4_w[0] = shift_poly_1(S_temp3_w[0]) ^ {9'b0, idata[39]};
							S_temp4_w[1] = shift_poly_2(S_temp3_w[1]) ^ {9'b0, idata[39]};
							S_temp4_w[2] = shift_poly_3(S_temp3_w[2]) ^ {9'b0, idata[39]};
							S_temp4_w[3] = shift_poly_4(S_temp3_w[3]) ^ {9'b0, idata[39]};
							S_temp4_w[4] = shift_poly_5(S_temp3_w[4]) ^ {9'b0, idata[39]};
							S_temp4_w[5] = shift_poly_6(S_temp3_w[5]) ^ {9'b0, idata[39]};
							S_temp4_w[6] = shift_poly_7(S_temp3_w[6]) ^ {9'b0, idata[39]};
							S_temp4_w[7] = shift_poly_8(S_temp3_w[7]) ^ {9'b0, idata[39]};
							S_temp5_w[0] = shift_poly_1(S_temp4_w[0]) ^ {9'b0, idata[31]};
							S_temp5_w[1] = shift_poly_2(S_temp4_w[1]) ^ {9'b0, idata[31]};
							S_temp5_w[2] = shift_poly_3(S_temp4_w[2]) ^ {9'b0, idata[31]};
							S_temp5_w[3] = shift_poly_4(S_temp4_w[3]) ^ {9'b0, idata[31]};
							S_temp5_w[4] = shift_poly_5(S_temp4_w[4]) ^ {9'b0, idata[31]};
							S_temp5_w[5] = shift_poly_6(S_temp4_w[5]) ^ {9'b0, idata[31]};
							S_temp5_w[6] = shift_poly_7(S_temp4_w[6]) ^ {9'b0, idata[31]};
							S_temp5_w[7] = shift_poly_8(S_temp4_w[7]) ^ {9'b0, idata[31]};
							S_temp6_w[0] = shift_poly_1(S_temp5_w[0]) ^ {9'b0, idata[23]};
							S_temp6_w[1] = shift_poly_2(S_temp5_w[1]) ^ {9'b0, idata[23]};
							S_temp6_w[2] = shift_poly_3(S_temp5_w[2]) ^ {9'b0, idata[23]};
							S_temp6_w[3] = shift_poly_4(S_temp5_w[3]) ^ {9'b0, idata[23]};
							S_temp6_w[4] = shift_poly_5(S_temp5_w[4]) ^ {9'b0, idata[23]};
							S_temp6_w[5] = shift_poly_6(S_temp5_w[5]) ^ {9'b0, idata[23]};
							S_temp6_w[6] = shift_poly_7(S_temp5_w[6]) ^ {9'b0, idata[23]};
							S_temp6_w[7] = shift_poly_8(S_temp5_w[7]) ^ {9'b0, idata[23]};
							S_temp7_w[0] = shift_poly_1(S_temp6_w[0]) ^ {9'b0, idata[15]};
							S_temp7_w[1] = shift_poly_2(S_temp6_w[1]) ^ {9'b0, idata[15]};
							S_temp7_w[2] = shift_poly_3(S_temp6_w[2]) ^ {9'b0, idata[15]};
							S_temp7_w[3] = shift_poly_4(S_temp6_w[3]) ^ {9'b0, idata[15]};
							S_temp7_w[4] = shift_poly_5(S_temp6_w[4]) ^ {9'b0, idata[15]};
							S_temp7_w[5] = shift_poly_6(S_temp6_w[5]) ^ {9'b0, idata[15]};
							S_temp7_w[6] = shift_poly_7(S_temp6_w[6]) ^ {9'b0, idata[15]};
							S_temp7_w[7] = shift_poly_8(S_temp6_w[7]) ^ {9'b0, idata[15]};
							S_w[0] = shift_poly_1(S_temp7_w[0]) ^ {9'b0, idata[7]};
							S_w[1] = shift_poly_2(S_temp7_w[1]) ^ {9'b0, idata[7]};
							S_w[2] = shift_poly_3(S_temp7_w[2]) ^ {9'b0, idata[7]};
							S_w[3] = shift_poly_4(S_temp7_w[3]) ^ {9'b0, idata[7]};
							S_w[4] = shift_poly_5(S_temp7_w[4]) ^ {9'b0, idata[7]};
							S_w[5] = shift_poly_6(S_temp7_w[5]) ^ {9'b0, idata[7]};
							S_w[6] = shift_poly_7(S_temp7_w[6]) ^ {9'b0, idata[7]};
							S_w[7] = shift_poly_8(S_temp7_w[7]) ^ {9'b0, idata[7]};
						end
						else if (cnt_r >= 7) begin
							S_temp1_w[0] = shift_poly_1(S_r[0]) ^ {9'b0, idata[63]};
							S_temp1_w[1] = shift_poly_2(S_r[1]) ^ {9'b0, idata[63]};
							S_temp1_w[2] = shift_poly_3(S_r[2]) ^ {9'b0, idata[63]};
							S_temp1_w[3] = shift_poly_4(S_r[3]) ^ {9'b0, idata[63]};
							S_temp1_w[4] = shift_poly_5(S_r[4]) ^ {9'b0, idata[63]};
							S_temp1_w[5] = shift_poly_6(S_r[5]) ^ {9'b0, idata[63]};
							S_temp1_w[6] = shift_poly_7(S_r[6]) ^ {9'b0, idata[63]};
							S_temp1_w[7] = shift_poly_8(S_r[7]) ^ {9'b0, idata[63]};
							S_temp2_w[0] = shift_poly_1(S_temp1_w[0]) ^ {9'b0, idata[55]};
							S_temp2_w[1] = shift_poly_2(S_temp1_w[1]) ^ {9'b0, idata[55]};
							S_temp2_w[2] = shift_poly_3(S_temp1_w[2]) ^ {9'b0, idata[55]};
							S_temp2_w[3] = shift_poly_4(S_temp1_w[3]) ^ {9'b0, idata[55]};
							S_temp2_w[4] = shift_poly_5(S_temp1_w[4]) ^ {9'b0, idata[55]};
							S_temp2_w[5] = shift_poly_6(S_temp1_w[5]) ^ {9'b0, idata[55]};
							S_temp2_w[6] = shift_poly_7(S_temp1_w[6]) ^ {9'b0, idata[55]};
							S_temp2_w[7] = shift_poly_8(S_temp1_w[7]) ^ {9'b0, idata[55]};
							S_temp3_w[0] = shift_poly_1(S_temp2_w[0]) ^ {9'b0, idata[47]};
							S_temp3_w[1] = shift_poly_2(S_temp2_w[1]) ^ {9'b0, idata[47]};
							S_temp3_w[2] = shift_poly_3(S_temp2_w[2]) ^ {9'b0, idata[47]};
							S_temp3_w[3] = shift_poly_4(S_temp2_w[3]) ^ {9'b0, idata[47]};
							S_temp3_w[4] = shift_poly_5(S_temp2_w[4]) ^ {9'b0, idata[47]};
							S_temp3_w[5] = shift_poly_6(S_temp2_w[5]) ^ {9'b0, idata[47]};
							S_temp3_w[6] = shift_poly_7(S_temp2_w[6]) ^ {9'b0, idata[47]};
							S_temp3_w[7] = shift_poly_8(S_temp2_w[7]) ^ {9'b0, idata[47]};
							S_temp4_w[0] = shift_poly_1(S_temp3_w[0]) ^ {9'b0, idata[39]};
							S_temp4_w[1] = shift_poly_2(S_temp3_w[1]) ^ {9'b0, idata[39]};
							S_temp4_w[2] = shift_poly_3(S_temp3_w[2]) ^ {9'b0, idata[39]};
							S_temp4_w[3] = shift_poly_4(S_temp3_w[3]) ^ {9'b0, idata[39]};
							S_temp4_w[4] = shift_poly_5(S_temp3_w[4]) ^ {9'b0, idata[39]};
							S_temp4_w[5] = shift_poly_6(S_temp3_w[5]) ^ {9'b0, idata[39]};
							S_temp4_w[6] = shift_poly_7(S_temp3_w[6]) ^ {9'b0, idata[39]};
							S_temp4_w[7] = shift_poly_8(S_temp3_w[7]) ^ {9'b0, idata[39]};
							S_temp5_w[0] = shift_poly_1(S_temp4_w[0]) ^ {9'b0, idata[31]};
							S_temp5_w[1] = shift_poly_2(S_temp4_w[1]) ^ {9'b0, idata[31]};
							S_temp5_w[2] = shift_poly_3(S_temp4_w[2]) ^ {9'b0, idata[31]};
							S_temp5_w[3] = shift_poly_4(S_temp4_w[3]) ^ {9'b0, idata[31]};
							S_temp5_w[4] = shift_poly_5(S_temp4_w[4]) ^ {9'b0, idata[31]};
							S_temp5_w[5] = shift_poly_6(S_temp4_w[5]) ^ {9'b0, idata[31]};
							S_temp5_w[6] = shift_poly_7(S_temp4_w[6]) ^ {9'b0, idata[31]};
							S_temp5_w[7] = shift_poly_8(S_temp4_w[7]) ^ {9'b0, idata[31]};
							S_temp6_w[0] = shift_poly_1(S_temp5_w[0]) ^ {9'b0, idata[23]};
							S_temp6_w[1] = shift_poly_2(S_temp5_w[1]) ^ {9'b0, idata[23]};
							S_temp6_w[2] = shift_poly_3(S_temp5_w[2]) ^ {9'b0, idata[23]};
							S_temp6_w[3] = shift_poly_4(S_temp5_w[3]) ^ {9'b0, idata[23]};
							S_temp6_w[4] = shift_poly_5(S_temp5_w[4]) ^ {9'b0, idata[23]};
							S_temp6_w[5] = shift_poly_6(S_temp5_w[5]) ^ {9'b0, idata[23]};
							S_temp6_w[6] = shift_poly_7(S_temp5_w[6]) ^ {9'b0, idata[23]};
							S_temp6_w[7] = shift_poly_8(S_temp5_w[7]) ^ {9'b0, idata[23]};
							S_temp7_w[0] = shift_poly_1(S_temp6_w[0]) ^ {9'b0, idata[15]};
							S_temp7_w[1] = shift_poly_2(S_temp6_w[1]) ^ {9'b0, idata[15]};
							S_temp7_w[2] = shift_poly_3(S_temp6_w[2]) ^ {9'b0, idata[15]};
							S_temp7_w[3] = shift_poly_4(S_temp6_w[3]) ^ {9'b0, idata[15]};
							S_temp7_w[4] = shift_poly_5(S_temp6_w[4]) ^ {9'b0, idata[15]};
							S_temp7_w[5] = shift_poly_6(S_temp6_w[5]) ^ {9'b0, idata[15]};
							S_temp7_w[6] = shift_poly_7(S_temp6_w[6]) ^ {9'b0, idata[15]};
							S_temp7_w[7] = shift_poly_8(S_temp6_w[7]) ^ {9'b0, idata[15]};
							S_w[0] = shift_poly_1(S_temp7_w[0]) ^ {9'b0, idata[7]};
							S_w[1] = shift_poly_2(S_temp7_w[1]) ^ {9'b0, idata[7]};
							S_w[2] = shift_poly_3(S_temp7_w[2]) ^ {9'b0, idata[7]};
							S_w[3] = shift_poly_4(S_temp7_w[3]) ^ {9'b0, idata[7]};
							S_w[4] = shift_poly_5(S_temp7_w[4]) ^ {9'b0, idata[7]};
							S_w[5] = shift_poly_6(S_temp7_w[5]) ^ {9'b0, idata[7]};
							S_w[6] = shift_poly_7(S_temp7_w[6]) ^ {9'b0, idata[7]};
							S_w[7] = shift_poly_8(S_temp7_w[7]) ^ {9'b0, idata[7]};
						end
					end
					default: state_w = S_IDLE;
				endcase
				if (mode_r) begin
					if (cnt_r >= 7) begin
						data_w[cnt_r - 7] = idata[7:0];
						data_w[cnt_r - 6] = idata[15:8];
						data_w[cnt_r - 5] = idata[23:16];
						data_w[cnt_r - 4] = idata[31:24];
						data_w[cnt_r - 3] = idata[39:32];
						data_w[cnt_r - 2] = idata[47:40];
						data_w[cnt_r - 1] = idata[55:48];
						data_w[cnt_r]     = idata[63:56];
						min1_w = min1_out;
						min2_w = min2_out;
						index1_w = index1_out;
						index2_w = index2_out;
					end
				end
				if (cnt_r > 7) begin
					cnt_w = cnt_r - 8;
				end
				else if (cnt_r == 7) begin
					ready_w = 0;
					cnt_w = cnt_r - 1;
				end else begin
					if(!mode_r) begin
						case (code_r)
						1: begin
							if (S1_S4_0) begin
								cnt_w = 2;
								state_w = S_OUT;
								finish_w = 1;
								odata_w = 1023;
							end
							else begin
								cnt_w = 0;
								state_w = S_BER;
								delta_w[0] = 1;
								delta_rho_w[0] = 1;
								// for (i = 0; i < 4; i = i + 1) $display("S%d = %b", i+1, S_r[i]);
								for (i = 1; i < 5; i = i + 1) begin
									delta_w[i] = 0;
									delta_rho_w[i] = 0;
								end
								rho_w = -1;
								l_rho_w = 0;
								l_w = 0;
								d_rho_w = 1;
								d_w = S_r[0];
							end
						end
						2: begin
							if (S1_S4_0) begin
								cnt_w = 2;
								state_w = S_OUT;
								finish_w = 1;
								odata_w = 1023;
							end
							else begin
								cnt_w = 0;
								state_w = S_BER;
								delta_w[0] = 1;
								delta_rho_w[0] = 1;
								// for (i = 0; i < 4; i = i + 1) $display("S%d = %b", i+1, S_r[i]);
								for (i = 1; i < 5; i = i + 1) begin
									delta_w[i] = 0;
									delta_rho_w[i] = 0;
								end
								rho_w = -1;
								l_rho_w = 0;
								l_w = 0;
								d_rho_w = 1;
								d_w = S_r[0];
							end
						end
						3: begin
							if (S1_S8_0) begin
								cnt_w = 4;
								state_w = S_OUT;
								finish_w = 1;
								odata_w = 1023;
							end
							else begin
								cnt_w = 0;
								state_w = S_BER;
								delta_w[0] = 1;
								delta_rho_w[0] = 1;
								// for (i = 0; i < 4; i = i + 1) $display("S%d = %b", i+1, S_r[i]);
								for (i = 1; i < 5; i = i + 1) begin
									delta_w[i] = 0;
									delta_rho_w[i] = 0;
								end
								rho_w = -1;
								l_rho_w = 0;
								l_w = 0;
								d_rho_w = 1;
								d_w = S_r[0];
							end
						end
						endcase	
					end
					else begin
						cnt_w = 2;
						state_w = S_OUT;
						finish_w = 1;
						odata_w = 1023;
						// for(i = 63; i >= 0; i = i - 1) begin
							// $display("data%d = %b", i, data_r[i]);
						// end
					end
				end
			end
			S_BER: begin
				// For μ = 0 to 2t-1:
				// If d_μ = 0:
				// 	Δ^(μ+1)(X) = Δ^(μ)(X)
				// 	l_μ+1 = l_μ
				
				// If d_μ ≠ 0:
				// 	Find ρ < μ with d_ρ ≠ 0 and (ρ - l_ρ) maximized
					
				// 	Δ^(μ+1)(X) = d_ρ · Δ^(μ)(X) + d_μ · X^(μ-ρ) · Δ^(ρ)(X)
				// 	l_μ+1 = max(l_μ, l_ρ + μ - ρ)
				
				// Compute next discrepancy:
				// 	d_μ+1 = S_μ+2 + Δ_1^(μ+1) · S_μ+1 + ... + Δ_l_{μ+1}^(μ+1) · S_μ+2-l_{μ+1}
				if (code_r != 3) begin
					delta_w[0] = S_r[0];
					delta_w[1] = element_mul(S_r[0], S_r[0]);
					delta_w[2] = S_r[2] ^ element_mul(S_r[0], S_r[1]);
					state_w = (mode_r) ? S_CHI_SOFT1 : S_CHI_HARD;
					case (code_r)
						1: begin
							cnt_w = 8;
						end
						2: begin
							cnt_w = 32;
						end
						default: begin
							cnt_w = 7;
						end
					endcase
				end
				else if (cnt_r < 16) begin
					cnt_w = cnt_r + 1;
					if (!cnt_r[0]) begin
						if (d_r == 0) begin
							for (i = 0; i < 5; i = i + 1) delta_w[i] = delta_r[i];
							l_w = l_r;
						end
						else begin
							cnt_temp_w = cnt_r >> 1;
							for (i = 0; i < 5; i = i + 1) begin
								temp1_w[i] = element_mul(d_rho_r, delta_r[i]);
							end
							for (i = 0; i < 5; i = i + 1) begin
								if ($signed(i[3:0]) - $signed(cnt_temp_w) + $signed(rho_r) >= 0) temp2_w[i] = delta_rho_r[$signed(i[3:0]) - $signed(cnt_temp_w) + $signed(rho_r)];
								else temp2_w[i] = 11'b0;
							end
							for (i = 0; i < 5; i = i + 1) begin
								temp3_w[i] = element_mul(d_r, temp2_w[i]);
							end
							for (i = 0; i < 5; i = i + 1) begin
								delta_w[i] = temp1_w[i] ^ temp3_w[i];
							end
							l_w = ($signed(l_r) > $signed(cnt_temp_w) + $signed(l_rho_r) - $signed(rho_r)) ? $signed(l_r) : $signed(cnt_temp_w) + $signed(l_rho_r) - $signed(rho_r);
						end
						if (d_r != 0 && $signed(cnt_temp_w) - $signed(l_r) > $signed(rho_r) - $signed(l_rho_r)) begin
							rho_w = cnt_temp_w;
							l_rho_w = l_r;
							d_rho_w = d_r;
							for (i = 0; i < 5; i = i + 1) begin
								delta_rho_w[i] = delta_r[i];
							end
						end
					end
					else begin
						cnt_temp_w = (cnt_r >> 1) + 3'd1;
						d_w = compute_d(cnt_temp_w, l_r);  // mu + 1 and l_{mu + 1}
					end
				end
				else begin
					state_w = (mode_r) ? S_CHI_SOFT1 : S_CHI_HARD;
					cnt_w = 128;
				end
			end
			S_CHI_HARD: begin
				case (code_r)
					1: begin
						if (cnt_r == 8) begin

						end
						else begin
							temp_root_w[0] = stage1_buff_r[0][0] ^ stage1_buff_r[0][1] ^ stage1_buff_r[0][2];
							if (cnt_r == 7) begin

							end
							else if (temp_root_w[0] == 0) begin
								temp_root_cnt_w[0] = root_cnt_r + 1;
								root_w[3'd1 - root_cnt_r] = (cnt_r[6:0] << 3) + 7'd7;
								if ((cnt_r[6:0] << 3) + 7'd7 != index1_out || (cnt_r[6:0] << 3) + 7'd7 != index2_out) begin
									corr_w[0] = corr_r[0] ^ 
								end
							end
							else temp_root_cnt_w[0] = root_cnt_r;
							for (i = 1; i < 8; i = i + 1) begin
								temp_root_w[i] = stage1_buff_r[i][0] ^ stage1_buff_r[i][1] ^ stage1_buff_r[i][2];
								if (temp_root_w[i] == 0) begin
									temp_root_cnt_w[i] = temp_root_cnt_w[i - 1] + 1;
									root_w[3'd1 - temp_root_cnt_w[i - 1]] = (cnt_r[6:0] << 3) + 7'd7 - i[6:0];
								end
								else begin
									temp_root_cnt_w[i] = temp_root_cnt_w[i - 1];
								end
							end
						end
						root_cnt_w = temp_root_cnt_w[7];
						cnt_w = cnt_r - 1;
						for (i = 0; i < 8; i = i + 1) begin
							for (j = 0; j < 3; j = j + 1) begin
								stage1_buff_w[i][j] = delta_poly_w[i][j];
							end
						end
						delta_w[0] = delta_r[0];
						delta_w[1] = shift_poly_8(delta_r[1]);
						delta_w[2] = shift_poly_8(shift_poly_8(delta_r[2]));
					end 
					
					2: begin
						if (cnt_r == 32) begin

						end
						else begin
							temp_root_w[0] = stage1_buff_r[0][0] ^ stage1_buff_r[0][1] ^ stage1_buff_r[0][2];
							if (cnt_r == 31) begin

							end
							else if (temp_root_w[0] == 0) begin
								temp_root_cnt_w[0] = root_cnt_r + 1;
								root_w[3'd1 - root_cnt_r] = (cnt_r[6:0] << 3) + 7'd7;
							end
							else temp_root_cnt_w[0] = root_cnt_r;
							for (i = 1; i < 8; i = i + 1) begin
								temp_root_w[i] = stage1_buff_r[i][0] ^ stage1_buff_r[i][1] ^ stage1_buff_r[i][2];
								if (temp_root_w[i] == 0) begin
									temp_root_cnt_w[i] = temp_root_cnt_w[i - 1] + 1;
									root_w[3'd1 - temp_root_cnt_w[i - 1]] = (cnt_r[6:0] << 3) + 7'd7 - i[6:0];
								end
								else begin
									temp_root_cnt_w[i] = temp_root_cnt_w[i - 1];
								end
							end
						end
						root_cnt_w = temp_root_cnt_w[7];
						cnt_w = cnt_r - 1;
						for (i = 0; i < 8; i = i + 1) begin
							for (j = 0; j < 3; j = j + 1) begin
								stage1_buff_w[i][j] = delta_poly_w[i][j];
							end
						end
						delta_w[0] = delta_r[0];
						delta_w[1] = shift_poly_8(delta_r[1]);
						delta_w[2] = shift_poly_8(shift_poly_8(delta_r[2]));
					end
					3: begin
						if (cnt_r == 128) begin

						end
						else begin
							temp_root_w[0] = stage1_buff_r[0][0] ^ stage1_buff_r[0][1] ^ stage1_buff_r[0][2] ^ stage1_buff_r[0][3] ^ stage1_buff_r[0][4];
							if (cnt_r == 127) begin

							end
							else if (temp_root_w[0] == 0) begin
								temp_root_cnt_w[0] = root_cnt_r + 1;
								root_w[3'd3 - root_cnt_r] = (cnt_r[6:0] << 3) + 7'd7;
							end
							else temp_root_cnt_w[0] = root_cnt_r;
							for (i = 1; i < 8; i = i + 1) begin
								temp_root_w[i] = stage1_buff_r[i][0] ^ stage1_buff_r[i][1] ^ stage1_buff_r[i][2] ^ stage1_buff_r[i][3] ^ stage1_buff_r[i][4];
								if (temp_root_w[i] == 0) begin
									temp_root_cnt_w[i] = temp_root_cnt_w[i - 1] + 1;
									root_w[3'd3 - temp_root_cnt_w[i - 1]] = (cnt_r[6:0] << 3) + 7'd7 - i[6:0];
								end
								else begin
									temp_root_cnt_w[i] = temp_root_cnt_w[i - 1];
								end
							end
						end
						root_cnt_w = temp_root_cnt_w[7];
						cnt_w = cnt_r - 1;
						for (i = 0; i < 8; i = i + 1) begin
							for (j = 0; j < 5; j = j + 1) begin
								stage1_buff_w[i][j] = delta_poly_w[i][j];
							end
						end
						delta_w[0] = delta_r[0];
						delta_w[1] = shift_poly_8(delta_r[1]);
						delta_w[2] = shift_poly_8(shift_poly_8(delta_r[2]));
						delta_w[3] = shift_poly_8(shift_poly_8(shift_poly_8(delta_r[3])));
						delta_w[4] = shift_poly_8(shift_poly_8(shift_poly_8(shift_poly_8(delta_r[4]))));
					end
					default: begin
					end
				endcase
				if (cnt_r == 0) begin
					state_w = S_OUT;
					cnt_w = 1;
					odata_w = root_r[0];
					finish_w = 1;
				end
			end
			S_CHI_SOFT1: begin
				case (code_r)
					1: begin
						if (cnt_r == 8) begin

						end
						else begin
							temp_root_w[0] = stage1_buff_r[0][0] ^ stage1_buff_r[0][1] ^ stage1_buff_r[0][2];
							if (cnt_r == 7) begin

							end
							else if (temp_root_w[0] == 0) begin
								temp_root_cnt_w[0] = root_cnt_r + 1;
								root_w[3'd1 - root_cnt_r] = (cnt_r[6:0] << 3) + 7'd7;
							end
							else temp_root_cnt_w[0] = root_cnt_r;
							for (i = 1; i < 8; i = i + 1) begin
								temp_root_w[i] = stage1_buff_r[i][0] ^ stage1_buff_r[i][1] ^ stage1_buff_r[i][2];
								if (temp_root_w[i] == 0) begin
									temp_root_cnt_w[i] = temp_root_cnt_w[i - 1] + 1;
									root_w[3'd1 - temp_root_cnt_w[i - 1]] = (cnt_r[6:0] << 3) + 7'd7 - i[6:0];
								end
								else begin
									temp_root_cnt_w[i] = temp_root_cnt_w[i - 1];
								end
							end
						end
						root_cnt_w = temp_root_cnt_w[7];
						cnt_w = cnt_r - 1;
						for (i = 0; i < 8; i = i + 1) begin
							for (j = 0; j < 3; j = j + 1) begin
								stage1_buff_w[i][j] = delta_poly_w[i][j];
							end
						end
						delta_w[0] = delta_r[0];
						delta_w[1] = shift_poly_8(delta_r[1]);
						delta_w[2] = shift_poly_8(shift_poly_8(delta_r[2]));
					end 
					
					2: begin
						if (cnt_r == 32) begin

						end
						else begin
							temp_root_w[0] = stage1_buff_r[0][0] ^ stage1_buff_r[0][1] ^ stage1_buff_r[0][2];
							if (cnt_r == 31) begin

							end
							else if (temp_root_w[0] == 0) begin
								temp_root_cnt_w[0] = root_cnt_r + 1;
								root_w[3'd1 - root_cnt_r] = (cnt_r[6:0] << 3) + 7'd7;
							end
							else temp_root_cnt_w[0] = root_cnt_r;
							for (i = 1; i < 8; i = i + 1) begin
								temp_root_w[i] = stage1_buff_r[i][0] ^ stage1_buff_r[i][1] ^ stage1_buff_r[i][2];
								if (temp_root_w[i] == 0) begin
									temp_root_cnt_w[i] = temp_root_cnt_w[i - 1] + 1;
									root_w[3'd1 - temp_root_cnt_w[i - 1]] = (cnt_r[6:0] << 3) + 7'd7 - i[6:0];
								end
								else begin
									temp_root_cnt_w[i] = temp_root_cnt_w[i - 1];
								end
							end
						end
						root_cnt_w = temp_root_cnt_w[7];
						cnt_w = cnt_r - 1;
						for (i = 0; i < 8; i = i + 1) begin
							for (j = 0; j < 3; j = j + 1) begin
								stage1_buff_w[i][j] = delta_poly_w[i][j];
							end
						end
						delta_w[0] = delta_r[0];
						delta_w[1] = shift_poly_8(delta_r[1]);
						delta_w[2] = shift_poly_8(shift_poly_8(delta_r[2]));
					end
					3: begin
						if (cnt_r == 128) begin

						end
						else begin
							temp_root_w[0] = stage1_buff_r[0][0] ^ stage1_buff_r[0][1] ^ stage1_buff_r[0][2] ^ stage1_buff_r[0][3] ^ stage1_buff_r[0][4];
							if (cnt_r == 127) begin

							end
							else if (temp_root_w[0] == 0) begin
								temp_root_cnt_w[0] = root_cnt_r + 1;
								root_w[3'd3 - root_cnt_r] = (cnt_r[6:0] << 3) + 7'd7;
							end
							else temp_root_cnt_w[0] = root_cnt_r;
							for (i = 1; i < 8; i = i + 1) begin
								temp_root_w[i] = stage1_buff_r[i][0] ^ stage1_buff_r[i][1] ^ stage1_buff_r[i][2] ^ stage1_buff_r[i][3] ^ stage1_buff_r[i][4];
								if (temp_root_w[i] == 0) begin
									temp_root_cnt_w[i] = temp_root_cnt_w[i - 1] + 1;
									root_w[3'd3 - temp_root_cnt_w[i - 1]] = (cnt_r[6:0] << 3) + 7'd7 - i[6:0];
								end
								else begin
									temp_root_cnt_w[i] = temp_root_cnt_w[i - 1];
								end
							end
						end
						root_cnt_w = temp_root_cnt_w[7];
						cnt_w = cnt_r - 1;
						for (i = 0; i < 8; i = i + 1) begin
							for (j = 0; j < 5; j = j + 1) begin
								stage1_buff_w[i][j] = delta_poly_w[i][j];
							end
						end
						delta_w[0] = delta_r[0];
						delta_w[1] = shift_poly_8(delta_r[1]);
						delta_w[2] = shift_poly_8(shift_poly_8(delta_r[2]));
						delta_w[3] = shift_poly_8(shift_poly_8(shift_poly_8(delta_r[3])));
						delta_w[4] = shift_poly_8(shift_poly_8(shift_poly_8(shift_poly_8(delta_r[4]))));

					end
					default: begin
					end
				endcase
				if (cnt_r == 0) begin
					state_w = S_OUT;
					cnt_w = 1;
					odata_w = root_r[0];
					finish_w = 1;
				end
			end
			S_CHI_SOFT2: begin
				
			end
			S_OUT: begin
				// $display("The error correcting function has below coefficients:");
				// for (i = 0; i < 5; i = i + 1) begin
				// 	$display("%b * X^ %d", delta_r[i], i[3:0]);
				// end
				// state_w = S_IDLE;
				// $display("index1: %d, val : %b", index1_r, min1_r);
				// $display("index2: %d, val : %b", index2_r, min2_r);
				odata_w = root_r[cnt_r[3:0]];
				if ((cnt_r == 2 && code_r != 3) || (cnt_r == 4 && code_r == 3)) begin
					finish_w = 0;
					odata_w = 0;
					state_w = S_IDLE;
					for (i = 0; i < 1024; i = i + 1) begin
						data_w[i] = 0;
					end
					for (i = 0; i < 8; i = i + 1) begin
						S_w[i] = 0;
						// power_w[i] = 0;
						// power2_w[i] = 0;
						// power3_w[i] = 0;
						// power4_w[i] = 0;
						for (j = 0; j < 5; j = j + 1) begin
							stage1_buff_w[i][j] = 0;
						end
					end
					for (i = 0; i < 5; i = i + 1) begin
						delta_w[i] = 0;
					end
					for (i = 0; i < 4; i = i + 1) begin
						root_w[i] = 0;
					end
					for (i = 0; i < 5; i = i + 1) begin
						delta_w[i] = 0;
						delta_rho_w[i] = 0;
					end
					d_w = 0;
					d_rho_w = 0;
					l_w = 0;
					l_rho_w = 0;
					rho_w = 0;
					cnt_w = 0;
					ready_w = 0;
					mode_w = 0;
					code_w = 0;
					root_cnt_w = 0;
				end
				else cnt_w = cnt_r + 1;
			end
		endcase
	end

	always @(posedge clk) begin
		if (!rstn) begin
			state_r <= S_IDLE;
			for (i = 0; i < 1024; i = i + 1) begin
				data_r[i] <= 0;
			end
			for (i = 0; i < 8; i = i + 1) begin
				S_r[i] <= 0;
				// power_r[i] <= 0;
				// power2_r[i] <= 0;
				// power3_r[i] <= 0;
				// power4_r[i] <= 0;
				for (j = 0; j < 5; j = j + 1) begin
					stage1_buff_r[i][j] <= 0;
				end
			end
			for (i = 0; i < 4; i = i + 1) begin
				root_r[i] <= 0;
				corr_r[i] <= 0;
			end
			for (i = 0; i < 5; i = i + 1) begin
				delta_r[i] <= 0;
				delta_rho_r[i] <= 0;
			end
			d_r <= 0;
			d_rho_r <= 0;
			l_r <= 0;
			l_rho_r <= 0;
			rho_r <= 0;
			cnt_r <= 0;
			ready_r <= 0;
			mode_r <= 0;
			code_r <= 0;
			root_cnt_r <= 0;
			odata_r <= 0;
			finish_r <= 0;
			min1_r <= 0;
			min2_r <= 0;
			index1_r <= 0;
			index2_r <= 0;
		end
		else begin
			cnt_r <= cnt_w;
			state_r <= state_w;
			ready_r <= ready_w;
			mode_r <= mode_w;
			code_r <= code_w;
			root_cnt_r <= root_cnt_w;
			min1_r <= min1_w;
			min2_r <= min2_w;
			index1_r <= index1_w;
			index2_r <= index2_w;
			if (load_en) begin
				for (i = 0; i < 1024; i = i + 1) begin
					data_r[i] <= data_w[i];
				end
			end
			if (syn_low_en) begin
				for (i = 0; i < 4; i = i + 1) begin
					S_r[i] <= S_w[i];
				end
			end
			if (syn_high_en) begin
				for (i = 0; i < 8; i = i + 1) begin
					S_r[i] <= S_w[i];
				end
			end
			if (ber_en) begin
				for (i = 0; i < 5; i = i + 1) begin
					delta_r[i] <= delta_w[i];
					delta_rho_r[i] <= delta_rho_w[i];
				end
				d_r <= d_w;
				d_rho_r <= d_rho_w;
				l_r <= l_w;
				l_rho_r <= l_rho_w;
				rho_r <= rho_w;
			end
			for (i = 0; i < 8; i = i + 1) begin
				// power_r[i] <= power_w[i];
				// power2_r[i] <= power2_w[i];
				// power3_r[i] <= power3_w[i];
				// power4_r[i] <= power4_w[i];
				for (j = 0; j < 5; j = j + 1) begin
					stage1_buff_r[i][j] <= stage1_buff_w[i][j];
				end
			end
			for (i = 0; i < 4; i = i + 1) begin
				root_r[i] <= root_w[i];
				corr_r[i] <= corr_w[i];
			end
			odata_r <= odata_w;
			finish_r <= finish_w;
		end
	end

	function automatic [10:0] poly_reduce_6;
		input [10:0] i_poly;
		begin
			poly_reduce_6 = i_poly;
			poly_reduce_6 = i_poly;
			poly_reduce_6[6] = 1'b0;  // Always clear bit 6
			poly_reduce_6[1] = i_poly[1] ^ i_poly[6];  // Conditional XOR
			poly_reduce_6[0] = i_poly[0] ^ i_poly[6];  // Conditional XOR
		end
		
	endfunction

	function automatic [10:0] poly_reduce_8;
		input [10:0] i_poly;
		begin
			poly_reduce_8 = i_poly;
			poly_reduce_8[8] = 1'b0;  // Always clear bit 8
			poly_reduce_8[0] = i_poly[0] ^ i_poly[8];  // Conditional XOR
			poly_reduce_8[2] = i_poly[2] ^ i_poly[8];  // Conditional XOR
			poly_reduce_8[3] = i_poly[3] ^ i_poly[8];  // Conditional XOR
			poly_reduce_8[4] = i_poly[4] ^ i_poly[8];  // Conditional XOR
		end
		
	endfunction

	function automatic [10:0] poly_reduce_10;
		input [10:0] i_poly;
		begin
			poly_reduce_10 = i_poly;
			poly_reduce_10[10] = 1'b0;  // Always clear bit 10
			poly_reduce_10[0] = i_poly[0] ^ i_poly[10];  // Conditional XOR
			poly_reduce_10[3] = i_poly[3] ^ i_poly[10];  // Conditional XOR
		end
		
	endfunction

	function automatic [10:0] shift_poly_1;
		input [9:0] i_poly;
		begin
			case (code_r)
				1: shift_poly_1 = poly_reduce_6(i_poly << 1);
				2: shift_poly_1 = poly_reduce_8(i_poly << 1);
				3: shift_poly_1 = poly_reduce_10(i_poly << 1); 
				default: shift_poly_1 = poly_reduce_6(i_poly << 1);
			endcase			
		end
	endfunction

	function automatic [10:0] shift_poly_2;
		input [9:0] i_poly;
		begin
			case (code_r)
				1: shift_poly_2 = poly_reduce_6(poly_reduce_6(i_poly << 1) << 1);
				2: shift_poly_2 = poly_reduce_8(poly_reduce_8(i_poly << 1) << 1);
				3: shift_poly_2 = poly_reduce_10(poly_reduce_10(i_poly << 1) << 1); 
				default: shift_poly_2 = poly_reduce_6(poly_reduce_6(i_poly << 1) << 1);
			endcase			
		end
	endfunction

	function automatic [10:0] shift_poly_3;
		input [9:0] i_poly;
		begin
			case (code_r)
				1: shift_poly_3 = poly_reduce_6(poly_reduce_6(poly_reduce_6(i_poly << 1) << 1) << 1);
				2: shift_poly_3 = poly_reduce_8(poly_reduce_8(poly_reduce_8(i_poly << 1) << 1) << 1);
				3: shift_poly_3 = poly_reduce_10(poly_reduce_10(poly_reduce_10(i_poly << 1) << 1) << 1); 
				default: shift_poly_3 = poly_reduce_6(poly_reduce_6(poly_reduce_6(i_poly << 1) << 1) << 1);
			endcase
		end
	endfunction

	function automatic [10:0] shift_poly_4;
		input [9:0] i_poly;
		begin
			case (code_r)
				1: shift_poly_4 = poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(i_poly << 1) << 1) << 1) << 1);
				2: shift_poly_4 = poly_reduce_8(poly_reduce_8(poly_reduce_8(poly_reduce_8(i_poly << 1) << 1) << 1) << 1);
				3: shift_poly_4 = poly_reduce_10(poly_reduce_10(poly_reduce_10(poly_reduce_10(i_poly << 1) << 1) << 1) << 1); 
				default: shift_poly_4 = poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(i_poly << 1) << 1) << 1) << 1);
			endcase
		end
	endfunction

	function automatic [10:0] shift_poly_5;
		input [9:0] i_poly;
		begin
			case (code_r)
				1: shift_poly_5 = poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(i_poly << 1) << 1) << 1) << 1) << 1);
				2: shift_poly_5 = poly_reduce_8(poly_reduce_8(poly_reduce_8(poly_reduce_8(poly_reduce_8(i_poly << 1) << 1) << 1) << 1) << 1);
				3: shift_poly_5 = poly_reduce_10(poly_reduce_10(poly_reduce_10(poly_reduce_10(poly_reduce_10(i_poly << 1) << 1) << 1) << 1) << 1); 
				default: shift_poly_5 = poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(i_poly << 1) << 1) << 1) << 1) << 1);
			endcase
		end
	endfunction

	function automatic[10:0] shift_poly_6;
        input [9:0] i_poly;
        begin
            case (code_r)
                1: shift_poly_6 = poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(i_poly << 1) << 1) << 1) << 1) << 1) << 1);
                2: shift_poly_6 = poly_reduce_8(poly_reduce_8(poly_reduce_8(poly_reduce_8(poly_reduce_8(poly_reduce_8(i_poly << 1) << 1) << 1) << 1) << 1) << 1);
                3: shift_poly_6 = poly_reduce_10(poly_reduce_10(poly_reduce_10(poly_reduce_10(poly_reduce_10(poly_reduce_10(i_poly << 1) << 1) << 1) << 1) << 1) << 1);
                default: shift_poly_6 = poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(i_poly << 1) << 1) << 1) << 1) << 1) << 1);
            endcase			
        end
    endfunction

    function automatic [10:0] shift_poly_7;
        input [9:0] i_poly;
        begin
            case (code_r)
                1: shift_poly_7 = poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(i_poly << 1) << 1) << 1) << 1) << 1) << 1) << 1);
                2: shift_poly_7 = poly_reduce_8(poly_reduce_8(poly_reduce_8(poly_reduce_8(poly_reduce_8(poly_reduce_8(poly_reduce_8(i_poly << 1) << 1) << 1) << 1) << 1) << 1) << 1);
                3: shift_poly_7 = poly_reduce_10(poly_reduce_10(poly_reduce_10(poly_reduce_10(poly_reduce_10(poly_reduce_10(poly_reduce_10(i_poly << 1) << 1) << 1) << 1) << 1) << 1) << 1);
                default: shift_poly_7 = poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(i_poly << 1) << 1) << 1) << 1) << 1) << 1) << 1);
            endcase			
        end
    endfunction

    function automatic [10:0] shift_poly_8;
        input [9:0] i_poly;
        begin
            case (code_r)
                1: shift_poly_8 = poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(i_poly << 1) << 1) << 1) << 1) << 1) << 1) << 1) << 1);
                2: shift_poly_8 = poly_reduce_8(poly_reduce_8(poly_reduce_8(poly_reduce_8(poly_reduce_8(poly_reduce_8(poly_reduce_8(poly_reduce_8(i_poly << 1) << 1) << 1) << 1) << 1) << 1) << 1) << 1);
                3: shift_poly_8 = poly_reduce_10(poly_reduce_10(poly_reduce_10(poly_reduce_10(poly_reduce_10(poly_reduce_10(poly_reduce_10(poly_reduce_10(i_poly << 1) << 1) << 1) << 1) << 1) << 1) << 1) << 1);
                default: shift_poly_8 = poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(poly_reduce_6(i_poly << 1) << 1) << 1) << 1) << 1) << 1) << 1) << 1);
            endcase			
        end
    endfunction

	function automatic [10:0] element_mul;
		input [10:0] i_element1;
		input [10:0] i_element2;
		reg [10:0] shift;
		reg [10:0] temp;
		integer i;
		begin
			shift = i_element1;
			element_mul = 0;
			temp = 0;
			for (i = 0; i < 11; i = i + 1) begin
				if (i_element2[i]) element_mul = element_mul ^ shift;
				temp = shift << 1;
				case (code_r)
					1: shift = poly_reduce_6(temp);
					2: shift = poly_reduce_8(temp);
					3: shift = poly_reduce_10(temp);
					default: shift = temp;
				endcase
			end
		end
	endfunction

	function automatic [10:0] compute_d;
		input [3:0]  i_mu_plus_1; // mu + 1
		input [3:0]  i_lu_plus_1; // l_{mu + 1}
		reg [3:0] cnt;
		integer i;
		begin
			cnt = i_mu_plus_1 + 1;
			compute_d = 0;
			for (i = 0; i < 5; i = i + 1) begin
				if ($signed(i[3:0]) <= $signed(i_lu_plus_1)) begin
					// S_{μ+2-i} = S_r[μ+1-i] = S_r[i_mu_plus_1 - i]
					if ($signed(i_mu_plus_1) >= $signed(i[3:0])) begin
						compute_d = compute_d ^ element_mul(delta_w[i], S_r[cnt - 1]);
						cnt = cnt - 1;
					end
				end
			end
		end
		
	endfunction

endmodule


module find_min (
	input [63:0] abs_data,
	input [9:0] cnt,
	input [7:0] min1_in,
	input [7:0] min2_in,
	input [9:0] index1_in,
	input [9:0] index2_in,
	output reg [7:0] min1_out,
	output reg [7:0] min2_out,
	output reg [9:0] index1_out,
	output reg [9:0] index2_out
);
	wire [7:0] abs_val [0:7];
    assign abs_val[0] = abs_data[63:56];
    assign abs_val[1] = abs_data[55:48];
    assign abs_val[2] = abs_data[47:40];
    assign abs_val[3] = abs_data[39:32];
    assign abs_val[4] = abs_data[31:24];
    assign abs_val[5] = abs_data[23:16];
    assign abs_val[6] = abs_data[15:8];
    assign abs_val[7] = abs_data[7:0];

	integer i;

	reg [7:0] stage_min1 [0:8];
    reg [7:0] stage_min2 [0:8];
    reg [9:0] stage_idx1 [0:8];
    reg [9:0] stage_idx2 [0:8];

	always @(*) begin
		stage_min1[0] = min1_in;
        stage_min2[0] = min2_in;
        stage_idx1[0] = index1_in;
        stage_idx2[0] = index2_in;
		for (i = 0; i < 8; i = i + 1) begin
            if (abs_val[i] < stage_min1[i]) begin
                stage_min1[i+1] = abs_val[i];
                stage_idx1[i+1] = cnt - i;
                stage_min2[i+1] = stage_min1[i];
                stage_idx2[i+1] = stage_idx1[i];
            end
            else if (abs_val[i] < stage_min2[i]) begin
                stage_min1[i+1] = stage_min1[i];
                stage_idx1[i+1] = stage_idx1[i];
                stage_min2[i+1] = abs_val[i];
                stage_idx2[i+1] = cnt - i;
            end
            else begin
                stage_min1[i+1] = stage_min1[i];
                stage_idx1[i+1] = stage_idx1[i];
                stage_min2[i+1] = stage_min2[i];
                stage_idx2[i+1] = stage_idx2[i];
            end
		end
		min1_out = stage_min1[8];
        min2_out = stage_min2[8];
        index1_out = stage_idx1[8];
        index2_out = stage_idx2[8];
	end
endmodule