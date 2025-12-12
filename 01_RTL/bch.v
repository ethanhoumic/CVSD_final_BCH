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
	localparam S_CORR_2  = 2;
	localparam S_BER_HARD  = 3;
	localparam S_BER_SOFT1 = 4;
	localparam S_BER_SOFT2 = 5;
	localparam S_CHI_HARD  = 6;
	localparam S_CHI_SOFT1 = 7;
	localparam S_CHI_SOFT2 = 8;
	localparam S_CORR_1 = 9;
	localparam S_OUT_HARD_BUFF = 10;
	localparam S_OUT_HARD  = 11;
	localparam S_OUT_SOFT_BUFF = 13;
	localparam S_OUT_SOFT_BUFF2 = 14;
	localparam S_OUT_SOFT  = 15;

	reg ready_r, ready_w;

	reg [3:0] state_r, state_w;
	reg [9:0] cnt_r, cnt_w;

	// input data
	reg [6:0] data_r [0:1023], data_w [0:1023];
	reg mode_r, mode_w;
	reg [1:0] code_r, code_w;

	// syndrome calculation
	reg [10:0] S_r [0:2][0:7], S_w [0:2][0:7];
	wire S1_S2_0 = (S_r[0][0] == 0 && S_r[0][1] == 0);
	wire S3_S4_0 = (S_r[0][2] == 0 && S_r[0][3] == 0);
	wire S5_S6_0 = (S_r[0][4] == 0 && S_r[0][5] == 0);
	wire S7_S8_0 = (S_r[0][6] == 0 && S_r[0][7] == 0);
	wire S1_S4_0 = S1_S2_0 && S3_S4_0;
	wire S5_S8_0 = S5_S6_0 && S7_S8_0;
	wire S1_S8_0 = S1_S4_0 && S5_S8_0;
	reg [9:0] S_temp1_w [0:7], S_temp2_w [0:7], S_temp3_w [0:7], S_temp4_w [0:7], S_temp5_w [0:7], S_temp6_w [0:7], S_temp7_w [0:7];
	
	reg [9:0] alpha_w [0:1][0:7], alpha_r[0:1][0:7];

	// soft decision
	wire [6:0] abs_idata [0:7];
	assign abs_idata[7] = idata[7] ? (~idata[6:0] + 1) : idata[6:0];
	assign abs_idata[6] = idata[15] ? (~idata[14:8] + 1) : idata[14:8];
	assign abs_idata[5] = idata[23] ? (~idata[22:16] + 1) : idata[22:16];
	assign abs_idata[4] = idata[31] ? (~idata[30:24] + 1) : idata[30:24];
	assign abs_idata[3] = idata[39] ? (~idata[38:32] + 1) : idata[38:32];
	assign abs_idata[2] = idata[47] ? (~idata[46:40] + 1) : idata[46:40];
	assign abs_idata[1] = idata[55] ? (~idata[54:48] + 1) : idata[54:48];
	assign abs_idata[0] = idata[63] ? (~idata[62:56] + 1) : idata[62:56];
	wire [6:0] abs_data0 = (((code_r == 3) && (cnt_r == 1023)) || ((code_r == 2) && (cnt_r == 255)) || ((code_r == 1) && (cnt_r == 63))) ? 7'b1111111 : abs_idata[0];
	wire [6:0] min1_cand, min2_cand;
	reg [6:0] min1_r, min2_r, min1_w, min2_w;
	reg [9:0] index1_r, index2_r, index1_w, index2_w;
	wire [9:0] index1_out, index2_out;
	wire [9:0] index1_cand, index2_cand;
	reg [9:0] index1_temp_w, index1_temp_r, index2_temp_w, index2_temp_r;

	wire [9:0] min1_out = 8'b0;
	wire [9:0] min2_out = 8'b0;

	swiss pipe_3_min(
		.clk(clk), .rst_n(rstn), .abs_data0(abs_data0), .abs_data1(abs_idata[1]), .abs_data2(abs_idata[2]), .abs_data3(abs_idata[3]), .abs_data4(abs_idata[4]), .abs_data5(abs_idata[5]), .abs_data6(abs_idata[6]), .abs_data7(abs_idata[7]), .cnt(cnt_r),
    	.min_1(min1_cand), .min_2(min2_cand), .index_1(index1_cand),. index_2(index2_cand)
	);

	//  ber algorithm
	reg [10:0] delta_r [0:2][0:6], delta_w [0:2][0:6], delta_rho_r [0:2][0:6], delta_rho_w [0:2][0:6], temp1_w [0:2][0:6], temp2_w [0:2][0:6], temp3_w [0:2][0:6];  // first index: big parallel number
	reg [10:0] d_r [0:2], d_w [0:2], d_rho_r[0:2], d_rho_w[0:2]; 
	reg [3:0]  l_r [0:2], l_w [0:2], l_rho_r [0:2], l_rho_w [0:2];
	reg [3:0]  rho_r [0:2], rho_w [0:2];
	reg [3:0]  cnt_temp_w [0:2];
	reg [2:0]  power_r [0:3], power_w [0:3];
	wire [4:0] ber_cnt_max_w = ((code_r != 3) ? 5'd8 : 5'd16);

	// chien search
	reg [2:0]  root_cnt_r [0:3], root_cnt_w [0:3];
	reg [10:0] temp_root_r [0:2][0:7], temp_root_w [0:2][0:7];
	reg [2:0]  temp_root_cnt_w [0:2][0:7];
	reg [9:0]  root_r [0:3][0:3], root_w [0:3][0:3];
	reg [10:0] delta_poly_w [0:2][0:7][0:4];    // first index: parallel number, second index: degree
	wire [9:0] root_index_w [0:7];
	wire [7:0] chien_cnt_max_w = ((code_r == 1) ? 8 : (code_r == 2) ? 32 : 128);

	genvar k, m;
	generate
		for (k = 0; k < 8; k = k + 1) begin
			assign root_index_w[k] = ((cnt_r << 3) + 7) - k;
		end
	endgenerate

	// correlation
	reg [9:0] corr_r[0:3], corr_w[0:3];
	reg [6:0] corr_cand_r [0:3][0:3], corr_cand_w [0:3][0:3];
	reg [1:0] corr_sel_r, corr_sel_w;
	reg [9:0] temp_corr_w [0:3][0:7];

	// output 
	reg [2:0] root_compact_w [0:4];
	reg [1:0] flag_w [0:4];
	reg [9:0] odata_r, odata_w;
	reg [9:0] flipped_stack_r [0:1], flipped_stack_w [0:1];
	reg [1:0] flipped_stack_ptr_r, flipped_stack_ptr_w;
	reg index1_invalid_r [0:3], index1_invalid_w [0:3], index2_invalid_r [0:3], index2_invalid_w [0:3];
	reg finish_r, finish_w;
	assign odata = odata_r;
	assign finish = finish_r;

	integer i, j, l;

	// clock gating
	wire data_low_en = mode_r;
	wire data_mid_en = mode_r && (code_r != 3);
	wire data_high_en = mode_r && (code_r == 3);
	// wire alpha_en = (state_r == S_CHI_SOFT1 || state_r == S_BER_SOFT1);
	// wire syn_low_en  = (state_r == S_LOAD || state_r == S_CHI_HARD || state_r == S_CHI_SOFT1);
	// wire syn_high_en = syn_low_en && (code_r == 3);
	// wire corr_en = (state_r == S_CORR_1 || state_r == S_CORR_2 || state_r == S_CHI_SOFT1);

	assign ready = ready_r;

	always @(*) begin
		for (i = 0; i < 8; i = i + 1) begin
			for (j = 0; j < 3; j = j + 1) begin
				delta_poly_w[j][i][0] = delta_r[j][0];
			end
		end
		for (i = 1; i < 5; i = i + 1) begin
			for (j = 0; j < 3; j = j + 1) begin
				delta_poly_w[j][0][i] = delta_r[j][i];
			end
		end
		case (code_r)
			1: begin
				for (i = 0; i < 3; i = i + 1) begin
					delta_poly_w[i][1][1] = shift_poly_1_6(delta_r[i][1]);
					delta_poly_w[i][1][2] = shift_poly_2_6(delta_r[i][2]);
					delta_poly_w[i][1][3] = shift_poly_3_6(delta_r[i][3]);
					delta_poly_w[i][1][4] = shift_poly_4_6(delta_r[i][4]);
					delta_poly_w[i][2][1] = shift_poly_2_6(delta_r[i][1]);
					delta_poly_w[i][2][2] = shift_poly_4_6(delta_r[i][2]);
					delta_poly_w[i][2][3] = shift_poly_6_6(delta_r[i][3]);
					delta_poly_w[i][2][4] = shift_poly_8_6(delta_r[i][4]);
					delta_poly_w[i][3][1] = shift_poly_3_6(delta_r[i][1]);
					delta_poly_w[i][3][2] = shift_poly_6_6(delta_r[i][2]);
					delta_poly_w[i][3][3] = shift_poly_9_6(delta_r[i][3]);
					delta_poly_w[i][3][4] = shift_poly_12_6(delta_r[i][4]);
					delta_poly_w[i][4][1] = shift_poly_4_6(delta_r[i][1]);
					delta_poly_w[i][4][2] = shift_poly_8_6(delta_r[i][2]);
					delta_poly_w[i][4][3] = shift_poly_12_6(delta_r[i][3]);
					delta_poly_w[i][4][4] = shift_poly_16_6(delta_r[i][4]);
					delta_poly_w[i][5][1] = shift_poly_5_6(delta_r[i][1]);
					delta_poly_w[i][5][2] = shift_poly_10_6(delta_r[i][2]);
					delta_poly_w[i][5][3] = shift_poly_15_6(delta_r[i][3]);
					delta_poly_w[i][5][4] = shift_poly_20_6(delta_r[i][4]);
					delta_poly_w[i][6][1] = shift_poly_6_6(delta_r[i][1]);
					delta_poly_w[i][6][2] = shift_poly_12_6(delta_r[i][2]);
					delta_poly_w[i][6][3] = shift_poly_18_6(delta_r[i][3]);
					delta_poly_w[i][6][4] = shift_poly_24_6(delta_r[i][4]);
					delta_poly_w[i][7][1] = shift_poly_7_6(delta_r[i][1]);
					delta_poly_w[i][7][2] = shift_poly_14_6(delta_r[i][2]);
					delta_poly_w[i][7][3] = shift_poly_21_6(delta_r[i][3]);
					delta_poly_w[i][7][4] = shift_poly_28_6(delta_r[i][4]);
				end
			end
			2: begin
				for (i = 0; i < 3; i = i + 1) begin
					delta_poly_w[i][1][1] = shift_poly_1_8(delta_r[i][1]);
					delta_poly_w[i][1][2] = shift_poly_2_8(delta_r[i][2]);
					delta_poly_w[i][1][3] = shift_poly_3_8(delta_r[i][3]);
					delta_poly_w[i][1][4] = shift_poly_4_8(delta_r[i][4]);
					delta_poly_w[i][2][1] = shift_poly_2_8(delta_r[i][1]);
					delta_poly_w[i][2][2] = shift_poly_4_8(delta_r[i][2]);
					delta_poly_w[i][2][3] = shift_poly_6_8(delta_r[i][3]);
					delta_poly_w[i][2][4] = shift_poly_8_8(delta_r[i][4]);
					delta_poly_w[i][3][1] = shift_poly_3_8(delta_r[i][1]);
					delta_poly_w[i][3][2] = shift_poly_6_8(delta_r[i][2]);
					delta_poly_w[i][3][3] = shift_poly_9_8(delta_r[i][3]);
					delta_poly_w[i][3][4] = shift_poly_12_8(delta_r[i][4]);
					delta_poly_w[i][4][1] = shift_poly_4_8(delta_r[i][1]);
					delta_poly_w[i][4][2] = shift_poly_8_8(delta_r[i][2]);
					delta_poly_w[i][4][3] = shift_poly_12_8(delta_r[i][3]);
					delta_poly_w[i][4][4] = shift_poly_16_8(delta_r[i][4]);
					delta_poly_w[i][5][1] = shift_poly_5_8(delta_r[i][1]);
					delta_poly_w[i][5][2] = shift_poly_10_8(delta_r[i][2]);
					delta_poly_w[i][5][3] = shift_poly_15_8(delta_r[i][3]);
					delta_poly_w[i][5][4] = shift_poly_20_8(delta_r[i][4]);
					delta_poly_w[i][6][1] = shift_poly_6_8(delta_r[i][1]);
					delta_poly_w[i][6][2] = shift_poly_12_8(delta_r[i][2]);
					delta_poly_w[i][6][3] = shift_poly_18_8(delta_r[i][3]);
					delta_poly_w[i][6][4] = shift_poly_24_8(delta_r[i][4]);
					delta_poly_w[i][7][1] = shift_poly_7_8(delta_r[i][1]);
					delta_poly_w[i][7][2] = shift_poly_14_8(delta_r[i][2]);
					delta_poly_w[i][7][3] = shift_poly_21_8(delta_r[i][3]);
					delta_poly_w[i][7][4] = shift_poly_28_8(delta_r[i][4]);
				end
			end
			3: begin
				for (i = 0; i < 3; i = i + 1) begin
					delta_poly_w[i][1][1] = shift_poly_1_10(delta_r[i][1]);
					delta_poly_w[i][1][2] = shift_poly_2_10(delta_r[i][2]);
					delta_poly_w[i][1][3] = shift_poly_3_10(delta_r[i][3]);
					delta_poly_w[i][1][4] = shift_poly_4_10(delta_r[i][4]);
					delta_poly_w[i][2][1] = shift_poly_2_10(delta_r[i][1]);
					delta_poly_w[i][2][2] = shift_poly_4_10(delta_r[i][2]);
					delta_poly_w[i][2][3] = shift_poly_6_10(delta_r[i][3]);
					delta_poly_w[i][2][4] = shift_poly_8_10(delta_r[i][4]);
					delta_poly_w[i][3][1] = shift_poly_3_10(delta_r[i][1]);
					delta_poly_w[i][3][2] = shift_poly_6_10(delta_r[i][2]);
					delta_poly_w[i][3][3] = shift_poly_9_10(delta_r[i][3]);
					delta_poly_w[i][3][4] = shift_poly_12_10(delta_r[i][4]);
					delta_poly_w[i][4][1] = shift_poly_4_10(delta_r[i][1]);
					delta_poly_w[i][4][2] = shift_poly_8_10(delta_r[i][2]);
					delta_poly_w[i][4][3] = shift_poly_12_10(delta_r[i][3]);
					delta_poly_w[i][4][4] = shift_poly_16_10(delta_r[i][4]);
					delta_poly_w[i][5][1] = shift_poly_5_10(delta_r[i][1]);
					delta_poly_w[i][5][2] = shift_poly_10_10(delta_r[i][2]);
					delta_poly_w[i][5][3] = shift_poly_15_10(delta_r[i][3]);
					delta_poly_w[i][5][4] = shift_poly_20_10(delta_r[i][4]);
					delta_poly_w[i][6][1] = shift_poly_6_10(delta_r[i][1]);
					delta_poly_w[i][6][2] = shift_poly_12_10(delta_r[i][2]);
					delta_poly_w[i][6][3] = shift_poly_18_10(delta_r[i][3]);
					delta_poly_w[i][6][4] = shift_poly_24_10(delta_r[i][4]);
					delta_poly_w[i][7][1] = shift_poly_7_10(delta_r[i][1]);
					delta_poly_w[i][7][2] = shift_poly_14_10(delta_r[i][2]);
					delta_poly_w[i][7][3] = shift_poly_21_10(delta_r[i][3]);
					delta_poly_w[i][7][4] = shift_poly_28_10(delta_r[i][4]);
				end
			end
			default: begin
				for (i = 0; i < 3; i = i + 1) begin
					delta_poly_w[i][1][1] = shift_poly_1_10(delta_r[i][1]);
					delta_poly_w[i][1][2] = shift_poly_2_10(delta_r[i][2]);
					delta_poly_w[i][1][3] = shift_poly_3_10(delta_r[i][3]);
					delta_poly_w[i][1][4] = shift_poly_4_10(delta_r[i][4]);
					delta_poly_w[i][2][1] = shift_poly_2_10(delta_r[i][1]);
					delta_poly_w[i][2][2] = shift_poly_4_10(delta_r[i][2]);
					delta_poly_w[i][2][3] = shift_poly_6_10(delta_r[i][3]);
					delta_poly_w[i][2][4] = shift_poly_8_10(delta_r[i][4]);
					delta_poly_w[i][3][1] = shift_poly_3_10(delta_r[i][1]);
					delta_poly_w[i][3][2] = shift_poly_6_10(delta_r[i][2]);
					delta_poly_w[i][3][3] = shift_poly_9_10(delta_r[i][3]);
					delta_poly_w[i][3][4] = shift_poly_12_10(delta_r[i][4]);
					delta_poly_w[i][4][1] = shift_poly_4_10(delta_r[i][1]);
					delta_poly_w[i][4][2] = shift_poly_8_10(delta_r[i][2]);
					delta_poly_w[i][4][3] = shift_poly_12_10(delta_r[i][3]);
					delta_poly_w[i][4][4] = shift_poly_16_10(delta_r[i][4]);
					delta_poly_w[i][5][1] = shift_poly_5_10(delta_r[i][1]);
					delta_poly_w[i][5][2] = shift_poly_10_10(delta_r[i][2]);
					delta_poly_w[i][5][3] = shift_poly_15_10(delta_r[i][3]);
					delta_poly_w[i][5][4] = shift_poly_20_10(delta_r[i][4]);
					delta_poly_w[i][6][1] = shift_poly_6_10(delta_r[i][1]);
					delta_poly_w[i][6][2] = shift_poly_12_10(delta_r[i][2]);
					delta_poly_w[i][6][3] = shift_poly_18_10(delta_r[i][3]);
					delta_poly_w[i][6][4] = shift_poly_24_10(delta_r[i][4]);
					delta_poly_w[i][7][1] = shift_poly_7_10(delta_r[i][1]);
					delta_poly_w[i][7][2] = shift_poly_14_10(delta_r[i][2]);
					delta_poly_w[i][7][3] = shift_poly_21_10(delta_r[i][3]);
					delta_poly_w[i][7][4] = shift_poly_28_10(delta_r[i][4]);
				end
			end
		endcase
	end

	// data_w
	always @(*) begin
		for (i = 0; i < 1024; i = i + 1) begin
			data_w[i] = data_r[i];
		end
		if ((state_r == S_LOAD) && (mode_r) && (cnt_r >= 7)) begin
			data_w[0] = abs_idata[7];
			data_w[1] = abs_idata[6];
			data_w[2] = abs_idata[5];
			data_w[3] = abs_idata[4];
			data_w[4] = abs_idata[3];
			data_w[5] = abs_idata[2];
			data_w[6] = abs_idata[1];
			data_w[7] = abs_idata[0];
			for (i = 8; i < 1024; i = i + 1) begin
				data_w[i] = data_r[i - 8];
			end
		end 
		else begin
			if ((state_r == S_CHI_SOFT1) || (state_r == S_CHI_SOFT2)) begin
				case (code_r)
					1: begin
						if ((cnt_r > 0) && (cnt_r < 9)) begin
							for (i = 8; i < 64; i = i + 1) begin
								data_w[i] = data_r[i - 8];
							end
							for(i = 0; i < 8; i = i + 1) begin
								data_w[i] = data_r[i + 56];
							end
						end
					end
					2: begin
						if ((cnt_r > 0) && (cnt_r < 33)) begin
							for (i = 8; i < 256; i = i + 1) begin
								data_w[i] = data_r[i - 8];
							end
							for(i = 0; i < 8; i = i + 1) begin
								data_w[i] = data_r[i + 248];
							end
						end
					end
					3: begin
						if ((cnt_r > 0) && (cnt_r < 129)) begin
							for (i = 8; i < 1024; i = i + 1) begin
								data_w[i] = data_r[i - 8];
							end
							for(i = 0; i < 8; i = i + 1) begin
								data_w[i] = data_r[i + 1016];
							end
						end
					end
					default: begin end
				endcase
			end
		end
	end

	always @(*) begin
		state_w = state_r;
		cnt_w = cnt_r;
		ready_w = ready_r;
		code_w = code_r;
		mode_w = mode_r;
		finish_w = finish_r;;
		odata_w = odata_r;
		corr_sel_w = corr_sel_r;
		flipped_stack_ptr_w = flipped_stack_ptr_r;
		for (i = 0; i < 2; i = i + 1) flipped_stack_w[i] = flipped_stack_r[i];
		for (i = 0; i < 5; i = i + 1) begin
			flag_w[i] = 0;
			root_compact_w[i] = 0;
			for (j = 0; j < 3; j = j + 1) begin
				temp1_w[j][i] = 11'b0;
				temp2_w[j][i] = 11'b0;
				temp3_w[j][i] = 11'b0;
			end
		end
		for (i = 0; i < 8; i = i + 1) begin
			for (j = 0; j < 3; j = j + 1) begin
				temp_root_w[j][i] = temp_root_r[j][i];
				temp_root_cnt_w[j][i] = 0;
				S_w[j][i] = S_r[j][i];
			end
		end 
		for (i = 0; i < 4; i = i + 1) begin
			power_w[i] = power_r[i];
			corr_w[i] = corr_r[i];
			root_cnt_w[i] = root_cnt_r[i];
			index1_invalid_w[i] = index1_invalid_r[i];
			index2_invalid_w[i] = index2_invalid_r[i];
			for (j = 0; j < 4; j = j + 1) begin
				root_w[i][j] = root_r[i][j];
				corr_cand_w[i][j] = corr_cand_r[i][j];
			end
		end
		for (i = 0; i < 3; i = i + 1) begin
			d_w[i] = d_r[i];
			d_rho_w[i] = d_rho_r[i];
			l_w[i] = l_r[i];
			l_rho_w[i] = l_rho_r[i];
			rho_w[i] = rho_r[i];
			cnt_temp_w[i] = 0;
			for (j = 0; j < 7; j = j + 1) begin
				delta_w[i][j] = delta_r[i][j];
				delta_rho_w[i][j] = delta_rho_r[i][j];
				temp1_w[i][j] = 11'b0;
				temp2_w[i][j] = 11'b0;
				temp3_w[i][j] = 11'b0;
			end
			for (j = 0; j < 8; j = j + 1) begin
				temp_root_w[i][j] = 1;
				temp_corr_w[i][j] = 0;
			end
		end
		case (state_r) // synopsys parallel_case
			S_IDLE: begin
				finish_w = 0;
				odata_w = 0;
				if (set && !ready_r) begin
					case (code) // synopsys parallel_case
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


			end  
			S_LOAD: begin
				case (code_r) // synopsys parallel_case
					1: begin
						if (cnt_r == 63) begin
							S_w[0][0] = {9'b0, idata[7]} ^ shift_poly_1_6({9'b0, idata[15]}) ^ shift_poly_2_6({9'b0, idata[23]}) ^ shift_poly_3_6({9'b0, idata[31]}) ^ shift_poly_4_6({9'b0, idata[39]}) ^ shift_poly_5_6({9'b0, idata[47]}) ^ shift_poly_6_6({9'b0, idata[55]});
							S_w[0][1] = {9'b0, idata[7]} ^ shift_poly_2_6({9'b0, idata[15]}) ^ shift_poly_4_6({9'b0, idata[23]}) ^ shift_poly_6_6({9'b0, idata[31]}) ^ shift_poly_8_6({9'b0, idata[39]}) ^ shift_poly_10_6({9'b0, idata[47]}) ^ shift_poly_12_6({9'b0, idata[55]});
							S_w[0][2] = {9'b0, idata[7]} ^ shift_poly_3_6({9'b0, idata[15]}) ^ shift_poly_6_6({9'b0, idata[23]}) ^ shift_poly_9_6({9'b0, idata[31]}) ^ shift_poly_12_6({9'b0, idata[39]}) ^ shift_poly_15_6({9'b0, idata[47]}) ^ shift_poly_18_6({9'b0, idata[55]});
							S_w[0][3] = {9'b0, idata[7]} ^ shift_poly_4_6({9'b0, idata[15]}) ^ shift_poly_8_6({9'b0, idata[23]}) ^ shift_poly_12_6({9'b0, idata[31]}) ^ shift_poly_16_6({9'b0, idata[39]}) ^ shift_poly_20_6({9'b0, idata[47]}) ^ shift_poly_24_6({9'b0, idata[55]});
						end
						else if (cnt_r >= 7) begin
							S_w[0][0] = {9'b0, idata[7]} ^ shift_poly_1_6({9'b0, idata[15]}) ^ shift_poly_2_6({9'b0, idata[23]}) ^ shift_poly_3_6({9'b0, idata[31]}) ^ shift_poly_4_6({9'b0, idata[39]}) ^ shift_poly_5_6({9'b0, idata[47]}) ^ shift_poly_6_6({9'b0, idata[55]}) ^ shift_poly_7_6({9'b0, idata[63]}) ^ shift_poly_8_6(S_r[0][0]);
							S_w[0][1] = {9'b0, idata[7]} ^ shift_poly_2_6({9'b0, idata[15]}) ^ shift_poly_4_6({9'b0, idata[23]}) ^ shift_poly_6_6({9'b0, idata[31]}) ^ shift_poly_8_6({9'b0, idata[39]}) ^ shift_poly_10_6({9'b0, idata[47]}) ^ shift_poly_12_6({9'b0, idata[55]}) ^ shift_poly_14_6({9'b0, idata[63]}) ^ shift_poly_16_6(S_r[0][1]);
							S_w[0][2] = {9'b0, idata[7]} ^ shift_poly_3_6({9'b0, idata[15]}) ^ shift_poly_6_6({9'b0, idata[23]}) ^ shift_poly_9_6({9'b0, idata[31]}) ^ shift_poly_12_6({9'b0, idata[39]}) ^ shift_poly_15_6({9'b0, idata[47]}) ^ shift_poly_18_6({9'b0, idata[55]}) ^ shift_poly_21_6({9'b0, idata[63]}) ^ shift_poly_24_6(S_r[0][2]);
							S_w[0][3] = {9'b0, idata[7]} ^ shift_poly_4_6({9'b0, idata[15]}) ^ shift_poly_8_6({9'b0, idata[23]}) ^ shift_poly_12_6({9'b0, idata[31]}) ^ shift_poly_16_6({9'b0, idata[39]}) ^ shift_poly_20_6({9'b0, idata[47]}) ^ shift_poly_24_6({9'b0, idata[55]}) ^ shift_poly_28_6({9'b0, idata[63]}) ^ shift_poly_32_6(S_r[0][3]);
						end
					end
					2: begin
						if (cnt_r == 255) begin
							S_w[0][0] = {9'b0, idata[7]} ^ shift_poly_1_8({9'b0, idata[15]}) ^ shift_poly_2_8({9'b0, idata[23]}) ^ shift_poly_3_8({9'b0, idata[31]}) ^ shift_poly_4_8({9'b0, idata[39]}) ^ shift_poly_5_8({9'b0, idata[47]}) ^ shift_poly_6_8({9'b0, idata[55]});
							S_w[0][1] = {9'b0, idata[7]} ^ shift_poly_2_8({9'b0, idata[15]}) ^ shift_poly_4_8({9'b0, idata[23]}) ^ shift_poly_6_8({9'b0, idata[31]}) ^ shift_poly_8_8({9'b0, idata[39]}) ^ shift_poly_10_8({9'b0, idata[47]}) ^ shift_poly_12_8({9'b0, idata[55]});
							S_w[0][2] = {9'b0, idata[7]} ^ shift_poly_3_8({9'b0, idata[15]}) ^ shift_poly_6_8({9'b0, idata[23]}) ^ shift_poly_9_8({9'b0, idata[31]}) ^ shift_poly_12_8({9'b0, idata[39]}) ^ shift_poly_15_8({9'b0, idata[47]}) ^ shift_poly_18_8({9'b0, idata[55]});
							S_w[0][3] = {9'b0, idata[7]} ^ shift_poly_4_8({9'b0, idata[15]}) ^ shift_poly_8_8({9'b0, idata[23]}) ^ shift_poly_12_8({9'b0, idata[31]}) ^ shift_poly_16_8({9'b0, idata[39]}) ^ shift_poly_20_8({9'b0, idata[47]}) ^ shift_poly_24_8({9'b0, idata[55]});
						end
						else if (cnt_r >= 7) begin
							S_w[0][0] = {9'b0, idata[7]} ^ shift_poly_1_8({9'b0, idata[15]}) ^ shift_poly_2_8({9'b0, idata[23]}) ^ shift_poly_3_8({9'b0, idata[31]}) ^ shift_poly_4_8({9'b0, idata[39]}) ^ shift_poly_5_8({9'b0, idata[47]}) ^ shift_poly_6_8({9'b0, idata[55]}) ^ shift_poly_7_8({9'b0, idata[63]}) ^ shift_poly_8_8(S_r[0][0]);
							S_w[0][1] = {9'b0, idata[7]} ^ shift_poly_2_8({9'b0, idata[15]}) ^ shift_poly_4_8({9'b0, idata[23]}) ^ shift_poly_6_8({9'b0, idata[31]}) ^ shift_poly_8_8({9'b0, idata[39]}) ^ shift_poly_10_8({9'b0, idata[47]}) ^ shift_poly_12_8({9'b0, idata[55]}) ^ shift_poly_14_8({9'b0, idata[63]}) ^ shift_poly_16_8(S_r[0][1]);
							S_w[0][2] = {9'b0, idata[7]} ^ shift_poly_3_8({9'b0, idata[15]}) ^ shift_poly_6_8({9'b0, idata[23]}) ^ shift_poly_9_8({9'b0, idata[31]}) ^ shift_poly_12_8({9'b0, idata[39]}) ^ shift_poly_15_8({9'b0, idata[47]}) ^ shift_poly_18_8({9'b0, idata[55]}) ^ shift_poly_21_8({9'b0, idata[63]}) ^ shift_poly_24_8(S_r[0][2]);
							S_w[0][3] = {9'b0, idata[7]} ^ shift_poly_4_8({9'b0, idata[15]}) ^ shift_poly_8_8({9'b0, idata[23]}) ^ shift_poly_12_8({9'b0, idata[31]}) ^ shift_poly_16_8({9'b0, idata[39]}) ^ shift_poly_20_8({9'b0, idata[47]}) ^ shift_poly_24_8({9'b0, idata[55]}) ^ shift_poly_28_8({9'b0, idata[63]}) ^ shift_poly_32_8(S_r[0][3]);
						end
					end
					3: begin
						if (cnt_r == 1023) begin
							S_w[0][0] = {9'b0, idata[7]} ^ shift_poly_1_10({9'b0, idata[15]}) ^ shift_poly_2_10({9'b0, idata[23]}) ^ shift_poly_3_10({9'b0, idata[31]}) ^ shift_poly_4_10({9'b0, idata[39]}) ^ shift_poly_5_10({9'b0, idata[47]}) ^ shift_poly_6_10({9'b0, idata[55]});
							S_w[0][1] = {9'b0, idata[7]} ^ shift_poly_2_10({9'b0, idata[15]}) ^ shift_poly_4_10({9'b0, idata[23]}) ^ shift_poly_6_10({9'b0, idata[31]}) ^ shift_poly_8_10({9'b0, idata[39]}) ^ shift_poly_10_10({9'b0, idata[47]}) ^ shift_poly_12_10({9'b0, idata[55]});
							S_w[0][2] = {9'b0, idata[7]} ^ shift_poly_3_10({9'b0, idata[15]}) ^ shift_poly_6_10({9'b0, idata[23]}) ^ shift_poly_9_10({9'b0, idata[31]}) ^ shift_poly_12_10({9'b0, idata[39]}) ^ shift_poly_15_10({9'b0, idata[47]}) ^ shift_poly_18_10({9'b0, idata[55]});
							S_w[0][3] = {9'b0, idata[7]} ^ shift_poly_4_10({9'b0, idata[15]}) ^ shift_poly_8_10({9'b0, idata[23]}) ^ shift_poly_12_10({9'b0, idata[31]}) ^ shift_poly_16_10({9'b0, idata[39]}) ^ shift_poly_20_10({9'b0, idata[47]}) ^ shift_poly_24_10({9'b0, idata[55]});
							S_w[0][4] = {9'b0, idata[7]} ^ shift_poly_5_10({9'b0, idata[15]}) ^ shift_poly_10_10({9'b0, idata[23]}) ^ shift_poly_15_10({9'b0, idata[31]}) ^ shift_poly_20_10({9'b0, idata[39]}) ^ shift_poly_25_10({9'b0, idata[47]}) ^ shift_poly_30_10({9'b0, idata[55]});
							S_w[0][5] = {9'b0, idata[7]} ^ shift_poly_6_10({9'b0, idata[15]}) ^ shift_poly_12_10({9'b0, idata[23]}) ^ shift_poly_18_10({9'b0, idata[31]}) ^ shift_poly_24_10({9'b0, idata[39]}) ^ shift_poly_30_10({9'b0, idata[47]}) ^ shift_poly_36_10({9'b0, idata[55]});
							S_w[0][6] = {9'b0, idata[7]} ^ shift_poly_7_10({9'b0, idata[15]}) ^ shift_poly_14_10({9'b0, idata[23]}) ^ shift_poly_21_10({9'b0, idata[31]}) ^ shift_poly_28_10({9'b0, idata[39]}) ^ shift_poly_35_10({9'b0, idata[47]}) ^ shift_poly_42_10({9'b0, idata[55]});
							S_w[0][7] = {9'b0, idata[7]} ^ shift_poly_8_10({9'b0, idata[15]}) ^ shift_poly_16_10({9'b0, idata[23]}) ^ shift_poly_24_10({9'b0, idata[31]}) ^ shift_poly_32_10({9'b0, idata[39]}) ^ shift_poly_40_10({9'b0, idata[47]}) ^ shift_poly_48_10({9'b0, idata[55]});
						end
						else if (cnt_r >= 7) begin
							S_w[0][0] = {9'b0, idata[7]} ^ shift_poly_1_10({9'b0, idata[15]}) ^ shift_poly_2_10({9'b0, idata[23]}) ^ shift_poly_3_10({9'b0, idata[31]}) ^ shift_poly_4_10({9'b0, idata[39]}) ^ shift_poly_5_10({9'b0, idata[47]}) ^ shift_poly_6_10({9'b0, idata[55]}) ^ shift_poly_7_10({9'b0, idata[63]}) ^ shift_poly_8_10(S_r[0][0]);
							S_w[0][1] = {9'b0, idata[7]} ^ shift_poly_2_10({9'b0, idata[15]}) ^ shift_poly_4_10({9'b0, idata[23]}) ^ shift_poly_6_10({9'b0, idata[31]}) ^ shift_poly_8_10({9'b0, idata[39]}) ^ shift_poly_10_10({9'b0, idata[47]}) ^ shift_poly_12_10({9'b0, idata[55]}) ^ shift_poly_14_10({9'b0, idata[63]}) ^ shift_poly_16_10(S_r[0][1]);
							S_w[0][2] = {9'b0, idata[7]} ^ shift_poly_3_10({9'b0, idata[15]}) ^ shift_poly_6_10({9'b0, idata[23]}) ^ shift_poly_9_10({9'b0, idata[31]}) ^ shift_poly_12_10({9'b0, idata[39]}) ^ shift_poly_15_10({9'b0, idata[47]}) ^ shift_poly_18_10({9'b0, idata[55]}) ^ shift_poly_21_10({9'b0, idata[63]}) ^ shift_poly_24_10(S_r[0][2]);
							S_w[0][3] = {9'b0, idata[7]} ^ shift_poly_4_10({9'b0, idata[15]}) ^ shift_poly_8_10({9'b0, idata[23]}) ^ shift_poly_12_10({9'b0, idata[31]}) ^ shift_poly_16_10({9'b0, idata[39]}) ^ shift_poly_20_10({9'b0, idata[47]}) ^ shift_poly_24_10({9'b0, idata[55]}) ^ shift_poly_28_10({9'b0, idata[63]}) ^ shift_poly_32_10(S_r[0][3]);
							S_w[0][4] = {9'b0, idata[7]} ^ shift_poly_5_10({9'b0, idata[15]}) ^ shift_poly_10_10({9'b0, idata[23]}) ^ shift_poly_15_10({9'b0, idata[31]}) ^ shift_poly_20_10({9'b0, idata[39]}) ^ shift_poly_25_10({9'b0, idata[47]}) ^ shift_poly_30_10({9'b0, idata[55]}) ^ shift_poly_35_10({9'b0, idata[63]}) ^ shift_poly_40_10(S_r[0][4]);
							S_w[0][5] = {9'b0, idata[7]} ^ shift_poly_6_10({9'b0, idata[15]}) ^ shift_poly_12_10({9'b0, idata[23]}) ^ shift_poly_18_10({9'b0, idata[31]}) ^ shift_poly_24_10({9'b0, idata[39]}) ^ shift_poly_30_10({9'b0, idata[47]}) ^ shift_poly_36_10({9'b0, idata[55]}) ^ shift_poly_42_10({9'b0, idata[63]}) ^ shift_poly_48_10(S_r[0][5]);
							S_w[0][6] = {9'b0, idata[7]} ^ shift_poly_7_10({9'b0, idata[15]}) ^ shift_poly_14_10({9'b0, idata[23]}) ^ shift_poly_21_10({9'b0, idata[31]}) ^ shift_poly_28_10({9'b0, idata[39]}) ^ shift_poly_35_10({9'b0, idata[47]}) ^ shift_poly_42_10({9'b0, idata[55]}) ^ shift_poly_49_10({9'b0, idata[63]}) ^ shift_poly_56_10(S_r[0][6]);
							S_w[0][7] = {9'b0, idata[7]} ^ shift_poly_8_10({9'b0, idata[15]}) ^ shift_poly_16_10({9'b0, idata[23]}) ^ shift_poly_24_10({9'b0, idata[31]}) ^ shift_poly_32_10({9'b0, idata[39]}) ^ shift_poly_40_10({9'b0, idata[47]}) ^ shift_poly_48_10({9'b0, idata[55]}) ^ shift_poly_56_10({9'b0, idata[63]}) ^ shift_poly_64_10(S_r[0][7]);
						end
					end
					default: state_w = S_IDLE;
				endcase
				if (cnt_r > 7) begin
					cnt_w = cnt_r - 8;
				end
				else if (cnt_r > 4) begin
					ready_w = 0;
					cnt_w = cnt_r - 1;
				end else begin
					if (S1_S4_0 && code_r != 3 || S1_S8_0 && code_r == 3) begin
						cnt_w = 1023;
						state_w = S_OUT_HARD;
						finish_w = 1;
						odata_w = 1023;
					end
					else begin
						cnt_w = 0;
						state_w = (mode_r) ? S_BER_SOFT1 : S_BER_HARD;
						delta_w[0][0] = 1;
						delta_rho_w[0][0] = 1;
						for (i = 1; i < 6; i = i + 1) begin
							delta_w[0][i] = 0;
							delta_rho_w[0][i] = 0;
						end
						rho_w[0] = -1;
						l_rho_w[0] = 0;
						l_w[0] = 0;
						d_rho_w[0] = 1;
						d_w[0] = S_r[0][0];
					end
				end
			end
			S_BER_HARD: begin
				if (code_r != 3) begin
					delta_w[0][0] = S_r[0][0];
					delta_w[0][1] = element_mul(S_r[0][0], S_r[0][0]);
					delta_w[0][2] = S_r[0][2] ^ element_mul(S_r[0][0], S_r[0][1]);
					state_w = (mode_r) ? S_CHI_SOFT1 : S_CHI_HARD;
					if (delta_w[0][2] != 0) power_w[0] = 2;
					else if (delta_w[0][1] != 0) power_w[0] = 1;
					else power_w[0] = 0;
					case (code_r) // synopsys parallel_case
						1: begin
							cnt_w = 8;
						end
						2: begin
							cnt_w = 32;
						end
						default: begin
							cnt_w = 8;
						end
					endcase
				end
				else if (cnt_r < 16) begin
					cnt_w = cnt_r + 1;
					if (!cnt_r[0]) begin
						if (d_r[0] == 0) begin
							for (i = 0; i < 5; i = i + 1) delta_w[0][i] = delta_r[0][i];
							l_w[0] = l_r[0];
						end
						else begin
							cnt_temp_w[0] = cnt_r >> 1;
							for (i = 0; i < 5; i = i + 1) temp1_w[0][i] = element_mul(d_rho_r[0], delta_r[0][i]);
							for (i = 0; i < 5; i = i + 1) begin
								if ($signed(i[3:0]) - $signed(cnt_temp_w[0]) + $signed(rho_r[0]) >= 0) temp2_w[0][i] = delta_rho_r[0][$signed(i[3:0]) - $signed(cnt_temp_w[0]) + $signed(rho_r[0])];
								else temp2_w[0][i] = 11'b0;
							end
							for (i = 0; i < 5; i = i + 1) temp3_w[0][i] = element_mul(d_r[0], temp2_w[0][i]);
							for (i = 0; i < 5; i = i + 1) delta_w[0][i] = temp1_w[0][i] ^ temp3_w[0][i];
							l_w[0] = ($signed(l_r[0]) > $signed(cnt_temp_w[0]) + $signed(l_rho_r[0]) - $signed(rho_r[0])) ? $signed(l_r[0]) : $signed(cnt_temp_w[0]) + $signed(l_rho_r[0]) - $signed(rho_r[0]);
						end
						if (d_r[0] != 0 && $signed(cnt_temp_w[0]) - $signed(l_r[0]) > $signed(rho_r[0]) - $signed(l_rho_r[0])) begin
							rho_w[0] = cnt_temp_w[0];
							l_rho_w[0] = l_r[0];
							d_rho_w[0] = d_r[0];
							for (i = 0; i < 5; i = i + 1) delta_rho_w[0][i] = delta_r[0][i];
						end
					end
					else begin
						cnt_temp_w[0] = (cnt_r >> 1) + 3'd1;
						d_w[0] = compute_d(cnt_temp_w[0], l_r[0], 2'd0);  // mu + 1 and l_{mu + 1}
					end
				end
				else begin
					if (delta_r[0][4] != 0) power_w[0] = 4;
					else if (delta_r[0][3] != 0) power_w[0] = 3;
					else if (delta_r[0][2] != 0) power_w[0] = 2;
					else if (delta_r[0][1] != 0) power_w[0] = 1;
					else power_w[0] = 0;
					state_w = (mode_r) ? S_CHI_SOFT1 : S_CHI_HARD;
					cnt_w = 128;
				end
			end
			S_BER_SOFT1: begin
				if (cnt_r < ber_cnt_max_w) begin
					cnt_w = cnt_r + 1;
					if (!cnt_r[0]) begin
						if (d_r[0] == 0) begin
							for (i = 0; i < 7; i = i + 1) delta_w[0][i] = delta_r[0][i];
							l_w[0] = l_r[0];
						end
						else begin
							cnt_temp_w[0] = cnt_r >> 1;
							for (i = 0; i < 7; i = i + 1) temp1_w[0][i] = element_mul(d_rho_r[0], delta_r[0][i]);
							for (i = 0; i < 7; i = i + 1) begin
								if ($signed(i[3:0]) - $signed(cnt_temp_w[0]) + $signed(rho_r[0]) >= 0) temp2_w[0][i] = delta_rho_r[0][$signed(i[3:0]) - $signed(cnt_temp_w[0]) + $signed(rho_r[0])];
								else temp2_w[0][i] = 11'b0;
							end
							for (i = 0; i < 7; i = i + 1) temp3_w[0][i] = element_mul(d_r[0], temp2_w[0][i]);
							for (i = 0; i < 7; i = i + 1) delta_w[0][i] = temp1_w[0][i] ^ temp3_w[0][i];
							l_w[0] = ($signed(l_r[0]) > $signed(cnt_temp_w[0]) + $signed(l_rho_r[0]) - $signed(rho_r[0])) ? $signed(l_r[0]) : $signed(cnt_temp_w[0]) + $signed(l_rho_r[0]) - $signed(rho_r[0]);
						end
						if (d_r[0] != 0 && $signed(cnt_temp_w[0]) - $signed(l_r[0]) > $signed(rho_r[0]) - $signed(l_rho_r[0])) begin
							rho_w[0] = cnt_temp_w[0];
							l_rho_w[0] = l_r[0];
							d_rho_w[0] = d_r[0];
							for (i = 0; i < 7; i = i + 1) delta_rho_w[0][i] = delta_r[0][i];
						end
					end
					else begin
						cnt_temp_w[0] = (cnt_r >> 1) + 3'd1;
						d_w[0] = compute_d(cnt_temp_w[0], l_r[0], 2'd0);  // mu + 1 and l_{mu + 1}
					end
				end
				else begin
					if (delta_r[0][6] != 0) power_w[0] = 6;
					else if (delta_r[0][5] != 0) power_w[0] = 5;
					else if (delta_r[0][4] != 0) power_w[0] = 4;
					else if (delta_r[0][3] != 0) power_w[0] = 3;
					else if (delta_r[0][2] != 0) power_w[0] = 2;
					else if (delta_r[0][1] != 0) power_w[0] = 1;
					else power_w[0] = 0;
					state_w = S_CHI_SOFT1;
					cnt_w = chien_cnt_max_w;
				end
			end
			S_BER_SOFT2: begin
				for (i = 0; i < 3; i = i + 1) begin
					if (cnt_r < ber_cnt_max_w) begin
						cnt_w = cnt_r + 1;
						if (!cnt_r[0]) begin
							if (d_r[i] == 0) begin
								for (j = 0; j < 7; j = j + 1) delta_w[i][j] = delta_r[i][j];
								l_w[i] = l_r[i];
							end
							else begin
								cnt_temp_w[i] = cnt_r >> 1;
								for (j = 0; j < 7; j = j + 1) temp1_w[i][j] = element_mul(d_rho_r[i], delta_r[i][j]);
								for (j = 0; j < 7; j = j + 1) begin
									if ($signed(j[3:0]) - $signed(cnt_temp_w[i]) + $signed(rho_r[i]) >= 0) temp2_w[i][j] = delta_rho_r[i][$signed(j[3:0]) - $signed(cnt_temp_w[i]) + $signed(rho_r[i])];
									else temp2_w[i][j] = 11'b0;
								end
								for (j = 0; j < 7; j = j + 1) temp3_w[i][j] = element_mul(d_r[i], temp2_w[i][j]);
								for (j = 0; j < 7; j = j + 1) delta_w[i][j] = temp1_w[i][j] ^ temp3_w[i][j];
								l_w[i] = ($signed(l_r[i]) > $signed(cnt_temp_w[i]) + $signed(l_rho_r[i]) - $signed(rho_r[i])) ? $signed(l_r[i]) : $signed(cnt_temp_w[i]) + $signed(l_rho_r[i]) - $signed(rho_r[i]);
							end
							if (d_r[i] != 0 && $signed(cnt_temp_w[i]) - $signed(l_r[i]) > $signed(rho_r[i]) - $signed(l_rho_r[i])) begin
								rho_w[i] = cnt_temp_w[i];
								l_rho_w[i] = l_r[i];
								d_rho_w[i] = d_r[i];
								for (j = 0; j < 7; j = j + 1) delta_rho_w[i][j] = delta_r[i][j];
							end
						end
						else begin
							cnt_temp_w[i] = (cnt_r >> 1) + 3'd1;
							d_w[i] = compute_d(cnt_temp_w[i], l_r[i], i[1:0]);  // mu + 1 and l_{mu + 1}
						end
					end
					else begin
						state_w = S_CHI_SOFT2;
						if (delta_r[i][6] != 0) power_w[i + 1] = 6;
						else if (delta_r[i][5] != 0) power_w[i + 1] = 5;
						else if (delta_r[i][4] != 0) power_w[i + 1] = 4;
						else if (delta_r[i][3] != 0) power_w[i + 1] = 3;
						else if (delta_r[i][2] != 0) power_w[i + 1] = 2;
						else if (delta_r[i][1] != 0) power_w[i + 1] = 1;
						else power_w[i + 1] = 0;
						cnt_w = chien_cnt_max_w;
					end
				end
			end
			S_CHI_HARD, S_CHI_SOFT1: begin
				for (i = 0; i < 8; i = i + 1) temp_root_w[0][i] = (code_r != 3) ? delta_poly_w[0][i][0] ^ delta_poly_w[0][i][1] ^ delta_poly_w[0][i][2] : delta_poly_w[0][i][0] ^ delta_poly_w[0][i][1] ^ delta_poly_w[0][i][2] ^ delta_poly_w[0][i][3] ^ delta_poly_w[0][i][4];
				if (cnt_r == chien_cnt_max_w) begin
					
				end
				else begin
					if (cnt_r == chien_cnt_max_w - 1) begin end
					else if (temp_root_r[0][0] == 0) begin
						root_w[0][root_cnt_r[0]] = root_index_w[0];
						temp_root_cnt_w[0][0] = root_cnt_r[0] + 1;
						corr_cand_w[0][root_cnt_r[0]] = data_r[7];
					end
					else temp_root_cnt_w[0][0] = root_cnt_r[0];
					for (i = 1; i < 8; i = i + 1) begin
						if (temp_root_r[0][i] == 0) begin
							temp_root_cnt_w[0][i] = temp_root_cnt_w[0][i - 1] + 1;
							root_w[0][temp_root_cnt_w[0][i - 1]] = root_index_w[i];
							corr_cand_w[0][temp_root_cnt_w[0][i - 1]] = data_r[7 - i];
						end
						else temp_root_cnt_w[0][i] = temp_root_cnt_w[0][i - 1];
					end
				end
				root_cnt_w[0] = temp_root_cnt_w[0][7];
				
				if (root_cnt_w[0] == 2 && code_r != 3 && !mode_r) begin
					state_w = S_OUT_HARD_BUFF;
					cnt_w = 1;
				end
				else if (root_cnt_w[0] == 4 && code_r == 3 && !mode_r) begin
					state_w = S_OUT_HARD_BUFF;
					cnt_w = 3;
				end
				else begin
					cnt_w = cnt_r - 1;
					case (code_r)
						1: begin
							delta_w[0][0] = delta_r[0][0];
							delta_w[0][1] = shift_poly_8_6(delta_r[0][1]);
							delta_w[0][2] = shift_poly_16_6(delta_r[0][2]);
						end 
						2: begin
							delta_w[0][0] = delta_r[0][0];
							delta_w[0][1] = shift_poly_8_8(delta_r[0][1]);
							delta_w[0][2] = shift_poly_16_8(delta_r[0][2]);
						end
						3: begin
							delta_w[0][0] = delta_r[0][0];
							delta_w[0][1] = shift_poly_8_10(delta_r[0][1]);
							delta_w[0][2] = shift_poly_16_10(delta_r[0][2]);
							delta_w[0][3] = shift_poly_24_10(delta_r[0][3]);
							delta_w[0][4] = shift_poly_32_10(delta_r[0][4]);
						end
						default: begin end
					endcase
				end
				if (cnt_r == 0) begin
					if (!mode_r) begin
						state_w = S_OUT_HARD_BUFF;
						cnt_w = root_cnt_w[0] - 1;
					end
					else begin
						state_w = S_BER_SOFT2;
						cnt_w = 0;
						for(i = 0; i < 8; i = i + 1) begin
							S_w[0][i] = S_r[0][i] ^ alpha_r[0][i];
							S_w[1][i] = S_r[0][i] ^ alpha_r[1][i];
							S_w[2][i] = S_w[0][i] ^ alpha_r[1][i];
						end
						for (i = 0; i < 3; i = i + 1) begin
							delta_w[i][0] = 1;
							delta_rho_w[i][0] = 1;
							rho_w[i] = -1;
							l_rho_w[i] = 0;
							l_w[i] = 0;
							d_rho_w[i] = 1;
							d_w[i] = S_w[i][0];
							for (j = 1; j < 6; j = j + 1) begin
								delta_w[i][j] = 0;
								delta_rho_w[i][j] = 0;
							end
						end
					end
				end
			end
			S_CHI_SOFT2: begin
				for (l = 0; l < 3; l = l + 1) begin
					for (i = 0; i < 8; i = i + 1) temp_root_w[l][i] = (code_r != 3) ? delta_poly_w[l][i][0] ^ delta_poly_w[l][i][1] ^ delta_poly_w[l][i][2] : delta_poly_w[l][i][0] ^ delta_poly_w[l][i][1] ^ delta_poly_w[l][i][2] ^ delta_poly_w[l][i][3] ^ delta_poly_w[l][i][4];
					if (cnt_r == chien_cnt_max_w) begin

					end
					else begin
						if (cnt_r == chien_cnt_max_w - 1) begin

						end
						else if (temp_root_r[l][0] == 0) begin
							temp_root_cnt_w[l][0] = root_cnt_r[l + 1] + 1;
							root_w[l + 1][root_cnt_r[l + 1]] = root_index_w[0];
							corr_cand_w[l + 1][root_cnt_r[l + 1]] = data_r[7];
						end
						else begin
							temp_root_cnt_w[l][0] = root_cnt_r[l + 1];
						end 
						for (i = 1; i < 8; i = i + 1) begin
							if (temp_root_r[l][i] == 0) begin
								temp_root_cnt_w[l][i] = temp_root_cnt_w[l][i - 1] + 1;
								root_w[l + 1][temp_root_cnt_w[l][i - 1]] = root_index_w[i];
								corr_cand_w[l + 1][temp_root_cnt_w[l][i - 1]] = data_r[7 - i];
							end
							else begin
								temp_root_cnt_w[l][i] = temp_root_cnt_w[l][i - 1];
							end
						end
					end
					root_cnt_w[l + 1] = temp_root_cnt_w[l][7];
					cnt_w = cnt_r - 1;
					case (code_r)
						1: begin
							delta_w[l][0] = delta_r[l][0];
							delta_w[l][1] = shift_poly_8_6(delta_r[l][1]);
							delta_w[l][2] = shift_poly_16_6(delta_r[l][2]);
						end
						2: begin
							delta_w[l][0] = delta_r[l][0];
							delta_w[l][1] = shift_poly_8_8(delta_r[l][1]);
							delta_w[l][2] = shift_poly_16_8(delta_r[l][2]);
						end
						3: begin
							delta_w[l][0] = delta_r[l][0];
							delta_w[l][1] = shift_poly_8_10(delta_r[l][1]);
							delta_w[l][2] = shift_poly_16_10(delta_r[l][2]);
							delta_w[l][3] = shift_poly_24_10(delta_r[l][3]);
							delta_w[l][4] = shift_poly_32_10(delta_r[l][4]);
						end
						default: begin end
					endcase
					
				end 
				if (cnt_r == 0) begin
					state_w = S_CORR_1;
					cnt_w = 0;
					corr_w[0] = 0;
					corr_w[1] = min1_r;
					corr_w[2] = min2_r;
					corr_w[3] = min1_r + min2_r;
				end
			end
			S_OUT_HARD_BUFF: begin
				odata_w = root_r[0][cnt_r[1:0]];
				finish_w = 1;
				state_w = S_OUT_HARD;
				cnt_w = cnt_r - 1;
			end
			S_OUT_HARD: begin
				odata_w = root_r[0][cnt_r[1:0]];
				if (cnt_r == 1023) begin
					finish_w = 0;
					odata_w = 0;
					state_w = S_IDLE;
					for (i = 0; i < 4; i = i + 1) begin
						power_w[0] = 0;
						root_w[0][i] = 0;
					end
					d_w[0] = 0;
					d_rho_w[0] = 0;
					l_w[0] = 0;
					l_rho_w[0] = 0;
					rho_w[0] = 0;
					root_cnt_w[0] = 0;
					cnt_w = 0;
					ready_w = 0;
					mode_w = 0;
					code_w = 0;
				end
				else cnt_w = cnt_r - 1;
			end
			S_CORR_1: begin
				if (cnt_r < root_cnt_r[0]) corr_w[0] = corr_r[0] + corr_cand_r[0][cnt_r[2:0]];
				else corr_w[0] = corr_r[0];

				if (cnt_r < root_cnt_r[1]) begin
					if (root_r[1][cnt_r[2:0]] == index1_r) corr_w[1] = corr_r[1] - min1_r;
					else corr_w[1] = corr_r[1] + corr_cand_r[1][cnt_r[2:0]];
				end
				else corr_w[1] = corr_r[1];

				if (cnt_r < root_cnt_r[2]) begin
					if (root_r[2][cnt_r[2:0]] == index2_r) corr_w[2] = corr_r[2] - min2_r;
					else corr_w[2] = corr_r[2] + corr_cand_r[2][cnt_r[2:0]];
				end
				else corr_w[2] = corr_r[2];

				if (cnt_r < root_cnt_r[3]) begin
					if (root_r[3][cnt_r[2:0]] == index1_r) corr_w[3] = corr_r[3] - min1_r;
					else if (root_r[3][cnt_r[2:0]] == index2_r) corr_w[3] = corr_r[3] - min2_r;
					else corr_w[3] = corr_r[3] + corr_cand_r[3][cnt_r[2:0]];
				end
				else corr_w[3] = corr_r[3];

				if ((cnt_r == 1 && code_r != 3) || (cnt_r == 3 || code_r == 3)) state_w = S_CORR_2;
				else begin
					state_w = state_r;
					cnt_w = cnt_r + 1;
				end
			end
			S_CORR_2: begin
				for (i = 0; i < 4; i = i + 1) begin
					if (root_cnt_r[i] < power_r[i] || (code_r != 3 && power_r[i] > 2) || (code_r == 3 && power_r[i] > 4) ) begin // decode failure
						corr_w[i] = 10'd1023;
					end
					else begin
						if (power_r[i] == 0) begin
							case (i[1:0])
								0: corr_w[0] = 0;
								1: corr_w[1] = min1_r;
								2: corr_w[2] = min2_r;
								3: corr_w[3] = min1_r + min2_r;
							endcase
						end
					end
				end
				corr_sel_w = compare_corr(corr_w[0], corr_w[1], corr_w[2], corr_w[3]);
				state_w = S_OUT_SOFT_BUFF;
			end
			S_OUT_SOFT_BUFF: begin
				flag_w[0] = 0;
				case (corr_sel_r)
					0: begin
						
					end 
					1: begin
						root_compact_w[0] = root_cnt_r[1];
						for (i = 1; i < 5; i = i + 1) begin
							if (i <= root_cnt_r[1]) begin
								if (root_r[1][i - 1] != index1_r) begin
									root_w[1][flag_w[i - 1]] = root_r[1][i - 1];
									flag_w[i] = flag_w[i - 1] + 1;
									root_compact_w[i] = root_compact_w[i - 1];
								end
								else begin
									index1_invalid_w[1] = 1;
									flag_w[i] = flag_w[i - 1];
									root_compact_w[i] = root_compact_w[i - 1] - 1;
								end
							end
							else root_compact_w[i] = root_compact_w[i - 1];
						end
						root_cnt_w[1] = root_compact_w[4];
					end
					2: begin
						root_compact_w[0] = root_cnt_r[2];
						for (i = 1; i < 5; i = i + 1) begin
							if (i <= root_cnt_r[2]) begin
								if (root_r[2][i - 1] != index2_r) begin
									root_w[2][flag_w[i - 1]] = root_r[2][i - 1];
									flag_w[i] = flag_w[i - 1] + 1;
									root_compact_w[i] = root_compact_w[i - 1];
								end
								else begin
									index2_invalid_w[2] = 1;
									flag_w[i] = flag_w[i - 1];
									root_compact_w[i] = root_compact_w[i - 1] - 1;
								end
							end
							else root_compact_w[i] = root_compact_w[i - 1];
						end
						root_cnt_w[2] = root_compact_w[4];
					end
					3: begin
						root_compact_w[0] = root_cnt_r[3];
						for (i = 1; i < 5; i = i + 1) begin
							if (i <= root_cnt_r[3]) begin
								if (root_r[3][i - 1] != index1_r && root_r[3][i - 1] != index2_r) begin
									root_w[3][flag_w[i - 1]] = root_r[3][i - 1];
									flag_w[i] = flag_w[i - 1] + 1;
									root_compact_w[i] = root_compact_w[i - 1];
								end
								else begin
									if (root_r[3][i - 1] == index1_r) index1_invalid_w[3] = 1;
									else index2_invalid_w[3] = 1;
									flag_w[i] = flag_w[i - 1];
									root_compact_w[i] = root_compact_w[i - 1] - 1;
								end
								end
							else root_compact_w[i] = root_compact_w[i - 1];
						end
						root_cnt_w[3] = root_compact_w[4];
					end
					default: begin end
				endcase
				state_w = S_OUT_SOFT_BUFF2;
			end
			S_OUT_SOFT_BUFF2: begin
				if (corr_sel_r == 0) begin
					if (power_r[0] == 0) cnt_w = 1023;
					else cnt_w = root_cnt_r[0] - 1;
					flipped_stack_ptr_w = 2;
				end 
				else if (corr_sel_r == 1)  begin
					if (power_r[1] == 0) cnt_w = 1023;
					else cnt_w = root_cnt_r[1] - 1;
					if (index1_invalid_r[1]) begin
						flipped_stack_ptr_w = 2;
					end
					else begin
						flipped_stack_w[1] = index1_r;
						flipped_stack_ptr_w = 1;
					end
				end 
				else if (corr_sel_r == 2) begin
					if (power_r[2] == 0) cnt_w = 1023;
					else cnt_w = root_cnt_r[2] - 1;
					if (index2_invalid_r[2]) begin
						flipped_stack_ptr_w = 2;
					end
					else begin
						flipped_stack_w[1] = index2_r;
						flipped_stack_ptr_w = 1;
					end
				end 
				else begin
					if (power_r[3] == 0) cnt_w = 1023;
					else cnt_w = root_cnt_r[3] - 1;
					if (index1_invalid_r[3] && index2_invalid_r[3]) begin
						flipped_stack_ptr_w = 2;
					end
					else if (index1_invalid_r[3] && !index2_invalid_r[3]) begin
						flipped_stack_w[1] = index2_r;
						flipped_stack_ptr_w = 1;
					end
					else if (!index1_invalid_r[3] && index2_invalid_r[3]) begin
						flipped_stack_w[1] = index1_r;
						flipped_stack_ptr_w = 1;
					end
					else begin
						if (index1_r < index2_r) begin
							flipped_stack_w[0] = index1_r;
							flipped_stack_w[1] = index2_r;
						end
						else begin
							flipped_stack_w[0] = index2_r;
							flipped_stack_w[1] = index1_r;
						end
						flipped_stack_ptr_w = 0;
					end
				end 
				state_w = S_OUT_SOFT;
			end
			S_OUT_SOFT: begin
				if (cnt_r != 1023) begin
					case (corr_sel_r) // synopsys parallel_case
						0: begin
							finish_w = 1;
							odata_w = root_r[0][cnt_r[1:0]];
							cnt_w = cnt_r - 1;
						end
						1: begin
							finish_w = 1;
							if (flipped_stack_ptr_r < 2) begin
								if (root_r[1][cnt_r[1:0]] < flipped_stack_r[flipped_stack_ptr_r]) begin
									odata_w = root_r[1][cnt_r[1:0]];
									cnt_w = cnt_r - 1;
								end
								else begin
									odata_w = flipped_stack_r[flipped_stack_ptr_r];
									flipped_stack_ptr_w = flipped_stack_ptr_r + 1;
								end
							end
							else begin
								odata_w = root_r[1][cnt_r[1:0]];
								cnt_w = cnt_r - 1;
							end
						end
						2: begin
							finish_w = 1;
							if (flipped_stack_ptr_r < 2) begin
								if (root_r[2][cnt_r[1:0]] < flipped_stack_r[flipped_stack_ptr_r]) begin
									odata_w = root_r[2][cnt_r[1:0]];
									cnt_w = cnt_r - 1;
								end
								else begin
									odata_w = flipped_stack_r[flipped_stack_ptr_r];
									flipped_stack_ptr_w = flipped_stack_ptr_r + 1;
								end
							end
							else begin
								odata_w = root_r[2][cnt_r[1:0]];
								cnt_w = cnt_r - 1;
							end
						end
						3: begin
							finish_w = 1;
							if (flipped_stack_ptr_r < 2) begin
								if (root_r[3][cnt_r[1:0]] < flipped_stack_r[flipped_stack_ptr_r]) begin
									odata_w = root_r[3][cnt_r[1:0]];
									cnt_w = cnt_r - 1;
								end
								else begin
									odata_w = flipped_stack_r[flipped_stack_ptr_r];
									flipped_stack_ptr_w = flipped_stack_ptr_r + 1;
								end
							end
							else begin
								odata_w = root_r[3][cnt_r[1:0]];
								cnt_w = cnt_r - 1;
							end
						end
						default: begin
							
						end
					endcase
				end
				else begin
					case (corr_sel_r) // synopsys parallel_case
						0: begin
							finish_w = 1;
						end
						1: begin
							finish_w = 1;
							if (flipped_stack_ptr_r < 2) begin
								odata_w = flipped_stack_r[flipped_stack_ptr_r];
								flipped_stack_ptr_w = flipped_stack_ptr_r + 1;
							end
						end
						2: begin
							finish_w = 1;
							if (flipped_stack_ptr_r < 2) begin
								odata_w = flipped_stack_r[flipped_stack_ptr_r];
								flipped_stack_ptr_w = flipped_stack_ptr_r + 1;
							end
						end
						3: begin
							finish_w = 1;
							if (flipped_stack_ptr_r < 2) begin
								odata_w = flipped_stack_r[flipped_stack_ptr_r];
								flipped_stack_ptr_w = flipped_stack_ptr_r + 1;
							end
						end
						default: begin
							
						end
					endcase
				end
				if (cnt_r == 1023 && flipped_stack_ptr_r == 2) begin
					finish_w = 0;
					odata_w = 0;
					state_w = S_IDLE;
					corr_sel_w = 0;
					flipped_stack_ptr_w = 0;
					flipped_stack_w[0] = 0;
					flipped_stack_w[1] = 0;
					for (i = 0; i < 4; i = i + 1) begin
						power_w[i] = 0;
						root_cnt_w[i] = 0;
						corr_w[i] = 0;
						index1_invalid_w[i] = 0;
						index2_invalid_w[i] = 0;
						for (j = 0; j < 4; j = j + 1) begin
							root_w[j][i] = 0;
						end
					end
					for (j = 0; j < 3; j = j + 1) begin
						d_w[j] = 0;
						d_rho_w[j] = 0;
						l_w[j] = 0;
						l_rho_w[j] = 0;
						rho_w[j] = 0;
					end
					flipped_stack_ptr_w = 0;
					cnt_w = 0;
					mode_w = 0;
					code_w = 0;
				end
			end
		endcase
	end

	// least reliability calc
	always @(*) begin
		min1_w = min1_r;
		min2_w = min2_r;
		index1_w = index1_r;
		index2_w = index2_r;
		if (state_r == S_LOAD) begin
			case (code_r)
				1: begin
					if (cnt_r == 47) begin
						min1_w = 7'b1111111;
						min2_w = 7'b1111111;
					end
					else if (cnt_r <= 39) begin	
						if (min1_cand < min1_r) begin
							min1_w = min1_cand;
							index1_w = index1_cand;
							min2_w = (min2_cand < min1_r) ? min2_cand : min1_r;
							index2_w = (min2_cand < min1_r) ? index2_cand : index1_r;
						end
						else if(min1_cand < min2_r) begin
							min2_w = min1_cand;
							index2_w = index1_cand;
						end
					end
				end
				2: begin
					if (cnt_r == 239) begin
						min1_w = 7'b1111111;
						min2_w = 7'b1111111;
					end
					else if (cnt_r <= 231) begin
						if (min1_cand < min1_r) begin
							min1_w = min1_cand;
							index1_w = index1_cand;
							min2_w = (min2_cand < min1_r) ? min2_cand : min1_r;
							index2_w = (min2_cand < min1_r) ? index2_cand : index1_r;
						end
						else if(min1_cand < min2_r) begin
							min2_w = min1_cand;
							index2_w = index1_cand;
						end
					end
				end
				3: begin
					if (cnt_r == 1007) begin
						min1_w = 7'b1111111;
						min2_w = 7'b1111111;
					end
					else if (cnt_r <= 999) begin
						if (min1_cand < min1_r) begin
							min1_w = min1_cand;
							index1_w = index1_cand;
							min2_w = (min2_cand < min1_r) ? min2_cand : min1_r;
							index2_w = (min2_cand < min1_r) ? index2_cand : index1_r;
						end
						else if(min1_cand < min2_r) begin
							min2_w = min1_cand;
							index2_w = index1_cand;
						end
					end
				end
			endcase
		end
	end

	// alpha calc
	always @(*) begin
		index1_temp_w = index1_temp_r;
		index2_temp_w = index2_temp_r;
		for (i = 0; i < 8; i = i + 1) begin
			alpha_w[0][i] = alpha_r[0][i];
			alpha_w[1][i] = alpha_r[1][i];
		end
		if (state_r == S_BER_SOFT1) begin
			index1_temp_w = index1_r;
			index2_temp_w = index2_r;
			for (i = 0; i < 8; i = i + 1) begin
				alpha_w[0][i] = 1;
				alpha_w[1][i] = 1;
			end	
		end 
		else if (state_r == S_CHI_SOFT1) begin
			index1_temp_w = (index1_temp_r > 7) ? index1_temp_r - 8 : 0;
			index2_temp_w = (index2_temp_r > 7) ? index2_temp_r - 8 : 0;
			case (code_r)
				1:begin
					if (index1_temp_r > 7) begin
						alpha_w[0][0] = shift_poly_8_6(alpha_r[0][0]);
						alpha_w[0][1] = shift_poly_16_6(alpha_r[0][1]);
						alpha_w[0][2] = shift_poly_24_6(alpha_r[0][2]);
						alpha_w[0][3] = shift_poly_32_6(alpha_r[0][3]);
					end
					else begin
						case (index1_temp_r[2:0])
							7: begin
								alpha_w[0][0] = shift_poly_7_6(alpha_r[0][0]);
								alpha_w[0][1] = shift_poly_14_6(alpha_r[0][1]);
								alpha_w[0][2] = shift_poly_21_6(alpha_r[0][2]);
								alpha_w[0][3] = shift_poly_28_6(alpha_r[0][3]);
							end
							6: begin
								alpha_w[0][0] = shift_poly_6_6(alpha_r[0][0]);
								alpha_w[0][1] = shift_poly_12_6(alpha_r[0][1]);
								alpha_w[0][2] = shift_poly_18_6(alpha_r[0][2]);
								alpha_w[0][3] = shift_poly_24_6(alpha_r[0][3]);
							end
							5: begin
								alpha_w[0][0] = shift_poly_5_6(alpha_r[0][0]);
								alpha_w[0][1] = shift_poly_10_6(alpha_r[0][1]);
								alpha_w[0][2] = shift_poly_15_6(alpha_r[0][2]);
								alpha_w[0][3] = shift_poly_20_6(alpha_r[0][3]);
							end
							4: begin
								alpha_w[0][0] = shift_poly_4_6(alpha_r[0][0]);
								alpha_w[0][1] = shift_poly_8_6(alpha_r[0][1]);
								alpha_w[0][2] = shift_poly_12_6(alpha_r[0][2]);
								alpha_w[0][3] = shift_poly_16_6(alpha_r[0][3]);
							end
							3: begin
								alpha_w[0][0] = shift_poly_3_6(alpha_r[0][0]);
								alpha_w[0][1] = shift_poly_6_6(alpha_r[0][1]);
								alpha_w[0][2] = shift_poly_9_6(alpha_r[0][2]);
								alpha_w[0][3] = shift_poly_12_6(alpha_r[0][3]);
							end
							2: begin
								alpha_w[0][0] = shift_poly_2_6(alpha_r[0][0]);
								alpha_w[0][1] = shift_poly_4_6(alpha_r[0][1]);
								alpha_w[0][2] = shift_poly_6_6(alpha_r[0][2]);
								alpha_w[0][3] = shift_poly_8_6(alpha_r[0][3]);
							end
							1: begin
								alpha_w[0][0] = shift_poly_1_6(alpha_r[0][0]);
								alpha_w[0][1] = shift_poly_2_6(alpha_r[0][1]);
								alpha_w[0][2] = shift_poly_3_6(alpha_r[0][2]);
								alpha_w[0][3] = shift_poly_4_6(alpha_r[0][3]);
							end
							default: begin
								alpha_w[0][0] = alpha_r[0][0];
								alpha_w[0][1] = alpha_r[0][1];
								alpha_w[0][2] = alpha_r[0][2];
								alpha_w[0][3] = alpha_r[0][3];
							end
						endcase
					end
					if (index2_temp_r > 7) begin
						alpha_w[1][0] = shift_poly_8_6(alpha_r[1][0]);
						alpha_w[1][1] = shift_poly_16_6(alpha_r[1][1]);
						alpha_w[1][2] = shift_poly_24_6(alpha_r[1][2]);
						alpha_w[1][3] = shift_poly_32_6(alpha_r[1][3]);
					end
					else begin
						case (index2_temp_r[2:0])
							7: begin
								alpha_w[1][0] = shift_poly_7_6(alpha_r[1][0]);
								alpha_w[1][1] = shift_poly_14_6(alpha_r[1][1]);
								alpha_w[1][2] = shift_poly_21_6(alpha_r[1][2]);
								alpha_w[1][3] = shift_poly_28_6(alpha_r[1][3]);
							end
							6: begin
								alpha_w[1][0] = shift_poly_6_6(alpha_r[1][0]);
								alpha_w[1][1] = shift_poly_12_6(alpha_r[1][1]);
								alpha_w[1][2] = shift_poly_18_6(alpha_r[1][2]);
								alpha_w[1][3] = shift_poly_24_6(alpha_r[1][3]);
							end
							5: begin
								alpha_w[1][0] = shift_poly_5_6(alpha_r[1][0]);
								alpha_w[1][1] = shift_poly_10_6(alpha_r[1][1]);
								alpha_w[1][2] = shift_poly_15_6(alpha_r[1][2]);
								alpha_w[1][3] = shift_poly_20_6(alpha_r[1][3]);
							end
							4: begin
								alpha_w[1][0] = shift_poly_4_6(alpha_r[1][0]);
								alpha_w[1][1] = shift_poly_8_6(alpha_r[1][1]);
								alpha_w[1][2] = shift_poly_12_6(alpha_r[1][2]);
								alpha_w[1][3] = shift_poly_16_6(alpha_r[1][3]);
							end
							3: begin
								alpha_w[1][0] = shift_poly_3_6(alpha_r[1][0]);
								alpha_w[1][1] = shift_poly_6_6(alpha_r[1][1]);
								alpha_w[1][2] = shift_poly_9_6(alpha_r[1][2]);
								alpha_w[1][3] = shift_poly_12_6(alpha_r[1][3]);
							end
							2: begin
								alpha_w[1][0] = shift_poly_2_6(alpha_r[1][0]);
								alpha_w[1][1] = shift_poly_4_6(alpha_r[1][1]);
								alpha_w[1][2] = shift_poly_6_6(alpha_r[1][2]);
								alpha_w[1][3] = shift_poly_8_6(alpha_r[1][3]);
							end
							1: begin
								alpha_w[1][0] = shift_poly_1_6(alpha_r[1][0]);
								alpha_w[1][1] = shift_poly_2_6(alpha_r[1][1]);
								alpha_w[1][2] = shift_poly_3_6(alpha_r[1][2]);
								alpha_w[1][3] = shift_poly_4_6(alpha_r[1][3]);
							end
							default: begin
								alpha_w[1][0] = alpha_r[1][0];
								alpha_w[1][1] = alpha_r[1][1];
								alpha_w[1][2] = alpha_r[1][2];
								alpha_w[1][3] = alpha_r[1][3];
							end
						endcase
					end
				end
				2: begin
					if (index1_temp_r > 7) begin
						alpha_w[0][0] = shift_poly_8_8(alpha_r[0][0]);
						alpha_w[0][1] = shift_poly_16_8(alpha_r[0][1]);
						alpha_w[0][2] = shift_poly_24_8(alpha_r[0][2]);
						alpha_w[0][3] = shift_poly_32_8(alpha_r[0][3]);
					end
					else begin
						case (index1_temp_r[2:0])
							7: begin
								alpha_w[0][0] = shift_poly_7_8(alpha_r[0][0]);
								alpha_w[0][1] = shift_poly_14_8(alpha_r[0][1]);
								alpha_w[0][2] = shift_poly_21_8(alpha_r[0][2]);
								alpha_w[0][3] = shift_poly_28_8(alpha_r[0][3]);
							end
							6: begin
								alpha_w[0][0] = shift_poly_6_8(alpha_r[0][0]);
								alpha_w[0][1] = shift_poly_12_8(alpha_r[0][1]);
								alpha_w[0][2] = shift_poly_18_8(alpha_r[0][2]);
								alpha_w[0][3] = shift_poly_24_8(alpha_r[0][3]);
							end
							5: begin
								alpha_w[0][0] = shift_poly_5_8(alpha_r[0][0]);
								alpha_w[0][1] = shift_poly_10_8(alpha_r[0][1]);
								alpha_w[0][2] = shift_poly_15_8(alpha_r[0][2]);
								alpha_w[0][3] = shift_poly_20_8(alpha_r[0][3]);
							end
							4: begin
								alpha_w[0][0] = shift_poly_4_8(alpha_r[0][0]);
								alpha_w[0][1] = shift_poly_8_8(alpha_r[0][1]);
								alpha_w[0][2] = shift_poly_12_8(alpha_r[0][2]);
								alpha_w[0][3] = shift_poly_16_8(alpha_r[0][3]);
							end
							3: begin
								alpha_w[0][0] = shift_poly_3_8(alpha_r[0][0]);
								alpha_w[0][1] = shift_poly_6_8(alpha_r[0][1]);
								alpha_w[0][2] = shift_poly_9_8(alpha_r[0][2]);
								alpha_w[0][3] = shift_poly_12_8(alpha_r[0][3]);
							end
							2: begin
								alpha_w[0][0] = shift_poly_2_8(alpha_r[0][0]);
								alpha_w[0][1] = shift_poly_4_8(alpha_r[0][1]);
								alpha_w[0][2] = shift_poly_6_8(alpha_r[0][2]);
								alpha_w[0][3] = shift_poly_8_8(alpha_r[0][3]);
							end
							1: begin
								alpha_w[0][0] = shift_poly_1_8(alpha_r[0][0]);
								alpha_w[0][1] = shift_poly_2_8(alpha_r[0][1]);
								alpha_w[0][2] = shift_poly_3_8(alpha_r[0][2]);
								alpha_w[0][3] = shift_poly_4_8(alpha_r[0][3]);
							end
							default: begin
								alpha_w[0][0] = alpha_r[0][0];
								alpha_w[0][1] = alpha_r[0][1];
								alpha_w[0][2] = alpha_r[0][2];
								alpha_w[0][3] = alpha_r[0][3];
							end
						endcase
					end
					if (index2_temp_r > 7) begin
						alpha_w[1][0] = shift_poly_8_8(alpha_r[1][0]);
						alpha_w[1][1] = shift_poly_16_8(alpha_r[1][1]);
						alpha_w[1][2] = shift_poly_24_8(alpha_r[1][2]);
						alpha_w[1][3] = shift_poly_32_8(alpha_r[1][3]);
					end
					else begin
						case (index2_temp_r[2:0])
							7: begin
								alpha_w[1][0] = shift_poly_7_8(alpha_r[1][0]);
								alpha_w[1][1] = shift_poly_14_8(alpha_r[1][1]);
								alpha_w[1][2] = shift_poly_21_8(alpha_r[1][2]);
								alpha_w[1][3] = shift_poly_28_8(alpha_r[1][3]);
							end
							6: begin
								alpha_w[1][0] = shift_poly_6_8(alpha_r[1][0]);
								alpha_w[1][1] = shift_poly_12_8(alpha_r[1][1]);
								alpha_w[1][2] = shift_poly_18_8(alpha_r[1][2]);
								alpha_w[1][3] = shift_poly_24_8(alpha_r[1][3]);
							end
							5: begin
								alpha_w[1][0] = shift_poly_5_8(alpha_r[1][0]);
								alpha_w[1][1] = shift_poly_10_8(alpha_r[1][1]);
								alpha_w[1][2] = shift_poly_15_8(alpha_r[1][2]);
								alpha_w[1][3] = shift_poly_20_8(alpha_r[1][3]);
							end
							4: begin
								alpha_w[1][0] = shift_poly_4_8(alpha_r[1][0]);
								alpha_w[1][1] = shift_poly_8_8(alpha_r[1][1]);
								alpha_w[1][2] = shift_poly_12_8(alpha_r[1][2]);
								alpha_w[1][3] = shift_poly_16_8(alpha_r[1][3]);
							end
							3: begin
								alpha_w[1][0] = shift_poly_3_8(alpha_r[1][0]);
								alpha_w[1][1] = shift_poly_6_8(alpha_r[1][1]);
								alpha_w[1][2] = shift_poly_9_8(alpha_r[1][2]);
								alpha_w[1][3] = shift_poly_12_8(alpha_r[1][3]);
							end
							2: begin
								alpha_w[1][0] = shift_poly_2_8(alpha_r[1][0]);
								alpha_w[1][1] = shift_poly_4_8(alpha_r[1][1]);
								alpha_w[1][2] = shift_poly_6_8(alpha_r[1][2]);
								alpha_w[1][3] = shift_poly_8_8(alpha_r[1][3]);
							end
							1: begin
								alpha_w[1][0] = shift_poly_1_8(alpha_r[1][0]);
								alpha_w[1][1] = shift_poly_2_8(alpha_r[1][1]);
								alpha_w[1][2] = shift_poly_3_8(alpha_r[1][2]);
								alpha_w[1][3] = shift_poly_4_8(alpha_r[1][3]);
							end
							default: begin
								alpha_w[1][0] = alpha_r[1][0];
								alpha_w[1][1] = alpha_r[1][1];
								alpha_w[1][2] = alpha_r[1][2];
								alpha_w[1][3] = alpha_r[1][3];
							end
						endcase
					end
				end
				3: begin
					if (index1_temp_r > 7) begin
						alpha_w[0][0] = shift_poly_8_10(alpha_r[0][0]);
						alpha_w[0][1] = shift_poly_16_10(alpha_r[0][1]);
						alpha_w[0][2] = shift_poly_24_10(alpha_r[0][2]);
						alpha_w[0][3] = shift_poly_32_10(alpha_r[0][3]);
						alpha_w[0][4] = shift_poly_40_10(alpha_r[0][4]);
						alpha_w[0][5] = shift_poly_48_10(alpha_r[0][5]);
						alpha_w[0][6] = shift_poly_56_10(alpha_r[0][6]);
						alpha_w[0][7] = shift_poly_64_10(alpha_r[0][7]);
					end
					else begin
						case (index1_temp_r[2:0])
							7: begin
								alpha_w[0][0] = shift_poly_7_10(alpha_r[0][0]);
								alpha_w[0][1] = shift_poly_14_10(alpha_r[0][1]);
								alpha_w[0][2] = shift_poly_21_10(alpha_r[0][2]);
								alpha_w[0][3] = shift_poly_28_10(alpha_r[0][3]);
								alpha_w[0][4] = shift_poly_35_10(alpha_r[0][4]);
								alpha_w[0][5] = shift_poly_42_10(alpha_r[0][5]);
								alpha_w[0][6] = shift_poly_49_10(alpha_r[0][6]);
								alpha_w[0][7] = shift_poly_56_10(alpha_r[0][7]);
							end
							6: begin
								alpha_w[0][0] = shift_poly_6_10(alpha_r[0][0]);
								alpha_w[0][1] = shift_poly_12_10(alpha_r[0][1]);
								alpha_w[0][2] = shift_poly_18_10(alpha_r[0][2]);
								alpha_w[0][3] = shift_poly_24_10(alpha_r[0][3]);
								alpha_w[0][4] = shift_poly_30_10(alpha_r[0][4]);
								alpha_w[0][5] = shift_poly_36_10(alpha_r[0][5]);
								alpha_w[0][6] = shift_poly_42_10(alpha_r[0][6]);
								alpha_w[0][7] = shift_poly_48_10(alpha_r[0][7]);
							end
							5: begin
								alpha_w[0][0] = shift_poly_5_10(alpha_r[0][0]);
								alpha_w[0][1] = shift_poly_10_10(alpha_r[0][1]);
								alpha_w[0][2] = shift_poly_15_10(alpha_r[0][2]);
								alpha_w[0][3] = shift_poly_20_10(alpha_r[0][3]);
								alpha_w[0][4] = shift_poly_25_10(alpha_r[0][4]);
								alpha_w[0][5] = shift_poly_30_10(alpha_r[0][5]);
								alpha_w[0][6] = shift_poly_35_10(alpha_r[0][6]);
								alpha_w[0][7] = shift_poly_40_10(alpha_r[0][7]);
							end
							4: begin
								alpha_w[0][0] = shift_poly_4_10(alpha_r[0][0]);
								alpha_w[0][1] = shift_poly_8_10(alpha_r[0][1]);
								alpha_w[0][2] = shift_poly_12_10(alpha_r[0][2]);
								alpha_w[0][3] = shift_poly_16_10(alpha_r[0][3]);
								alpha_w[0][4] = shift_poly_20_10(alpha_r[0][4]);
								alpha_w[0][5] = shift_poly_24_10(alpha_r[0][5]);
								alpha_w[0][6] = shift_poly_28_10(alpha_r[0][6]);
								alpha_w[0][7] = shift_poly_32_10(alpha_r[0][7]);
							end
							3: begin
								alpha_w[0][0] = shift_poly_3_10(alpha_r[0][0]);
								alpha_w[0][1] = shift_poly_6_10(alpha_r[0][1]);
								alpha_w[0][2] = shift_poly_9_10(alpha_r[0][2]);
								alpha_w[0][3] = shift_poly_12_10(alpha_r[0][3]);
								alpha_w[0][4] = shift_poly_15_10(alpha_r[0][4]);
								alpha_w[0][5] = shift_poly_18_10(alpha_r[0][5]);
								alpha_w[0][6] = shift_poly_21_10(alpha_r[0][6]);
								alpha_w[0][7] = shift_poly_24_10(alpha_r[0][7]);
							end
							2: begin
								alpha_w[0][0] = shift_poly_2_10(alpha_r[0][0]);
								alpha_w[0][1] = shift_poly_4_10(alpha_r[0][1]);
								alpha_w[0][2] = shift_poly_6_10(alpha_r[0][2]);
								alpha_w[0][3] = shift_poly_8_10(alpha_r[0][3]);
								alpha_w[0][4] = shift_poly_10_10(alpha_r[0][4]);
								alpha_w[0][5] = shift_poly_12_10(alpha_r[0][5]);
								alpha_w[0][6] = shift_poly_14_10(alpha_r[0][6]);
								alpha_w[0][7] = shift_poly_16_10(alpha_r[0][7]);
							end
							1: begin
								alpha_w[0][0] = shift_poly_1_10(alpha_r[0][0]);
								alpha_w[0][1] = shift_poly_2_10(alpha_r[0][1]);
								alpha_w[0][2] = shift_poly_3_10(alpha_r[0][2]);
								alpha_w[0][3] = shift_poly_4_10(alpha_r[0][3]);
								alpha_w[0][4] = shift_poly_5_10(alpha_r[0][4]);
								alpha_w[0][5] = shift_poly_6_10(alpha_r[0][5]);
								alpha_w[0][6] = shift_poly_7_10(alpha_r[0][6]);
								alpha_w[0][7] = shift_poly_8_10(alpha_r[0][7]);
							end
							default: begin
								alpha_w[0][0] = alpha_r[0][0];
								alpha_w[0][1] = alpha_r[0][1];
								alpha_w[0][2] = alpha_r[0][2];
								alpha_w[0][3] = alpha_r[0][3];
								alpha_w[0][4] = alpha_r[0][4];
								alpha_w[0][5] = alpha_r[0][5];
								alpha_w[0][6] = alpha_r[0][6];
								alpha_w[0][7] = alpha_r[0][7];
							end
						endcase
					end
					if (index2_temp_r > 7) begin
						alpha_w[1][0] = shift_poly_8_10(alpha_r[1][0]);
						alpha_w[1][1] = shift_poly_16_10(alpha_r[1][1]);
						alpha_w[1][2] = shift_poly_24_10(alpha_r[1][2]);
						alpha_w[1][3] = shift_poly_32_10(alpha_r[1][3]);
						alpha_w[1][4] = shift_poly_40_10(alpha_r[1][4]);
						alpha_w[1][5] = shift_poly_48_10(alpha_r[1][5]);
						alpha_w[1][6] = shift_poly_56_10(alpha_r[1][6]);
						alpha_w[1][7] = shift_poly_64_10(alpha_r[1][7]);
					end
					else begin
						case (index2_temp_r[2:0])
							7: begin
								alpha_w[1][0] = shift_poly_7_10(alpha_r[1][0]);
								alpha_w[1][1] = shift_poly_14_10(alpha_r[1][1]);
								alpha_w[1][2] = shift_poly_21_10(alpha_r[1][2]);
								alpha_w[1][3] = shift_poly_28_10(alpha_r[1][3]);
								alpha_w[1][4] = shift_poly_35_10(alpha_r[1][4]);
								alpha_w[1][5] = shift_poly_42_10(alpha_r[1][5]);
								alpha_w[1][6] = shift_poly_49_10(alpha_r[1][6]);
								alpha_w[1][7] = shift_poly_56_10(alpha_r[1][7]);
							end
							6: begin
								alpha_w[1][0] = shift_poly_6_10(alpha_r[1][0]);
								alpha_w[1][1] = shift_poly_12_10(alpha_r[1][1]);
								alpha_w[1][2] = shift_poly_18_10(alpha_r[1][2]);
								alpha_w[1][3] = shift_poly_24_10(alpha_r[1][3]);
								alpha_w[1][4] = shift_poly_30_10(alpha_r[1][4]);
								alpha_w[1][5] = shift_poly_36_10(alpha_r[1][5]);
								alpha_w[1][6] = shift_poly_42_10(alpha_r[1][6]);
								alpha_w[1][7] = shift_poly_48_10(alpha_r[1][7]);
							end
							5: begin
								alpha_w[1][0] = shift_poly_5_10(alpha_r[1][0]);
								alpha_w[1][1] = shift_poly_10_10(alpha_r[1][1]);
								alpha_w[1][2] = shift_poly_15_10(alpha_r[1][2]);
								alpha_w[1][3] = shift_poly_20_10(alpha_r[1][3]);
								alpha_w[1][4] = shift_poly_25_10(alpha_r[1][4]);
								alpha_w[1][5] = shift_poly_30_10(alpha_r[1][5]);
								alpha_w[1][6] = shift_poly_35_10(alpha_r[1][6]);
								alpha_w[1][7] = shift_poly_40_10(alpha_r[1][7]);
							end
							4: begin
								alpha_w[1][0] = shift_poly_4_10(alpha_r[1][0]);
								alpha_w[1][1] = shift_poly_8_10(alpha_r[1][1]);
								alpha_w[1][2] = shift_poly_12_10(alpha_r[1][2]);
								alpha_w[1][3] = shift_poly_16_10(alpha_r[1][3]);
								alpha_w[1][4] = shift_poly_20_10(alpha_r[1][4]);
								alpha_w[1][5] = shift_poly_24_10(alpha_r[1][5]);
								alpha_w[1][6] = shift_poly_28_10(alpha_r[1][6]);
								alpha_w[1][7] = shift_poly_32_10(alpha_r[1][7]);
							end
							3: begin
								alpha_w[1][0] = shift_poly_3_10(alpha_r[1][0]);
								alpha_w[1][1] = shift_poly_6_10(alpha_r[1][1]);
								alpha_w[1][2] = shift_poly_9_10(alpha_r[1][2]);
								alpha_w[1][3] = shift_poly_12_10(alpha_r[1][3]);
								alpha_w[1][4] = shift_poly_15_10(alpha_r[1][4]);
								alpha_w[1][5] = shift_poly_18_10(alpha_r[1][5]);
								alpha_w[1][6] = shift_poly_21_10(alpha_r[1][6]);
								alpha_w[1][7] = shift_poly_24_10(alpha_r[1][7]);
							end
							2: begin
								alpha_w[1][0] = shift_poly_2_10(alpha_r[1][0]);
								alpha_w[1][1] = shift_poly_4_10(alpha_r[1][1]);
								alpha_w[1][2] = shift_poly_6_10(alpha_r[1][2]);
								alpha_w[1][3] = shift_poly_8_10(alpha_r[1][3]);
								alpha_w[1][4] = shift_poly_10_10(alpha_r[1][4]);
								alpha_w[1][5] = shift_poly_12_10(alpha_r[1][5]);
								alpha_w[1][6] = shift_poly_14_10(alpha_r[1][6]);
								alpha_w[1][7] = shift_poly_16_10(alpha_r[1][7]);
							end
							1: begin
								alpha_w[1][0] = shift_poly_1_10(alpha_r[1][0]);
								alpha_w[1][1] = shift_poly_2_10(alpha_r[1][1]);
								alpha_w[1][2] = shift_poly_3_10(alpha_r[1][2]);
								alpha_w[1][3] = shift_poly_4_10(alpha_r[1][3]);
								alpha_w[1][4] = shift_poly_5_10(alpha_r[1][4]);
								alpha_w[1][5] = shift_poly_6_10(alpha_r[1][5]);
								alpha_w[1][6] = shift_poly_7_10(alpha_r[1][6]);
								alpha_w[1][7] = shift_poly_8_10(alpha_r[1][7]);
							end
							default: begin
								alpha_w[1][0] = alpha_r[1][0];
								alpha_w[1][1] = alpha_r[1][1];
								alpha_w[1][2] = alpha_r[1][2];
								alpha_w[1][3] = alpha_r[1][3];
								alpha_w[1][4] = alpha_r[1][4];
								alpha_w[1][5] = alpha_r[1][5];
								alpha_w[1][6] = alpha_r[1][6];
								alpha_w[1][7] = alpha_r[1][7];
							end
						endcase
					end
				end
			default: begin end
			endcase
		end
	end

	// gated
	wire gated_data_low_en = data_low_en | (!rstn);
	wire gated_data_mid_en = data_mid_en | (!rstn);
	wire gated_data_high_en = data_high_en | (!rstn);
	wire data_low_clk = gated_data_low_en & clk;
	wire data_mid_clk = gated_data_mid_en & clk;
	wire data_high_clk = gated_data_high_en & clk;

	always @(posedge data_low_clk) begin
		if (!rstn) begin
			for (i = 0; i < 64; i = i + 1) begin
				data_r[i] <= 0;
			end
		end
		else begin
			for (i = 0; i < 64; i = i + 1) begin
				data_r[i] <= data_w[i];
			end
		end
	end

	always @(posedge data_mid_clk) begin
		if (!rstn) begin
			for (i = 64; i < 256; i = i + 1) begin
				data_r[i] <= 0;
			end
		end
		else begin
			for (i = 64; i < 256; i = i + 1) begin
				data_r[i] <= data_w[i];
			end
		end
	end

	always @(posedge data_high_clk) begin
		if (!rstn) begin
			for (i = 256; i < 1024; i = i + 1) begin
				data_r[i] <= 0;
			end
		end
		else begin
			for (i = 256; i < 1024; i = i + 1) begin
				data_r[i] <= data_w[i];
			end
		end
	end


	always @(posedge clk) begin
		if (!rstn) begin
			state_r <= S_IDLE;
			for (i = 0; i < 8; i = i + 1)begin
				alpha_r[0][i] <= 0;
				alpha_r[1][i] <= 0;
			end
			for (i = 0; i < 4; i = i + 1) begin
				power_r[i] <= 0;
				corr_r[i] <= 0;
				root_cnt_r[i] <= 0;
				index1_invalid_r[i] <= 0;
				index2_invalid_r[i] <= 0;
				for (j = 0; j < 4; j = j + 1) begin
					root_r[i][j] <= 0;
					corr_cand_r[i][j] <= 0;
				end
			end
			cnt_r <= 0;
			ready_r <= 0;
			mode_r <= 0;
			code_r <= 0;
			odata_r <= 0;
			finish_r <= 0;
			min1_r <= 0;
			min2_r <= 0;
			index1_r <= 0;
			index2_r <= 0;
			index1_temp_r <= 0;
			index2_temp_r <= 0;
			corr_sel_r <= 0;
			flipped_stack_ptr_r <= 0;
			for (i = 0; i < 2; i = i + 1) flipped_stack_r[i] <= 0;
			for (i = 0; i < 3; i = i + 1) begin
				d_r[i] <= 0;
				d_rho_r[i] <= 0;
				l_r[i] <= 0;
				l_rho_r[i] <= 0;
				rho_r[i] <= 0;				
				for (j = 0; j < 7; j = j + 1) begin
					delta_r[i][j] <= 0;
					delta_rho_r[i][j] <= 0;
				end
				for (j = 0; j < 8; j = j + 1) begin
					temp_root_r[i][j] <= 0;
					S_r[i][j] <= 0;
				end
			end
		end
		else begin
			cnt_r <= cnt_w;
			state_r <= state_w;
			ready_r <= ready_w;
			mode_r <= mode_w;
			code_r <= code_w;
			min1_r <= min1_w;
			min2_r <= min2_w;
			index1_r <= index1_w;
			index2_r <= index2_w;
			index1_temp_r <= index1_temp_w;
			index2_temp_r <= index2_temp_w;
			corr_sel_r <= corr_sel_w;
			flipped_stack_ptr_r <= flipped_stack_ptr_w;
			for (i = 4; i < 8; i = i + 1) begin
				for (j = 0; j < 3; j = j + 1) begin
					S_r[j][i] <= S_w[j][i];
				end
			end
			for (i = 0; i < 4; i = i + 1) begin
				for (j = 0; j < 3; j = j + 1) begin
					S_r[j][i] <= S_w[j][i];
				end
			end
			for (i = 0; i < 3; i = i + 1) begin
				d_r[i] <= d_w[i];
				d_rho_r[i] <= d_rho_w[i];
				l_r[i] <= l_w[i];
				l_rho_r[i] <= l_rho_w[i];
				rho_r[i] <= rho_w[i];
				for (j = 0; j < 7; j = j + 1) begin
					delta_r[i][j] <= delta_w[i][j];
					delta_rho_r[i][j] <= delta_rho_w[i][j];
				end
			end
			for (i = 0; i < 2; i = i + 1) begin
				flipped_stack_r[i] <= flipped_stack_w[i];
			end
			for (i = 0; i < 8; i = i + 1)begin
				alpha_r[0][i] <= alpha_w[0][i];
				alpha_r[1][i] <= alpha_w[1][i];
			end
			for (i = 0; i < 4; i = i + 1) begin
				power_r[i] <= power_w[i];
				corr_r[i] <= corr_w[i];
				root_cnt_r[i] <= root_cnt_w[i];
				index1_invalid_r[i] <= index1_invalid_w[i];
				index2_invalid_r[i] <= index2_invalid_w[i];
				for (j = 0; j < 4; j = j + 1) begin
					root_r[i][j] <= root_w[i][j];
				end
			end
			for (i = 0; i < 4; i = i + 1) begin
				for (j = 0; j < 4; j = j + 1) begin
					corr_cand_r[i][j] <= corr_cand_w[i][j];
				end
			end
			for (i = 0; i < 3; i = i + 1) begin
				for (j = 0; j < 8; j = j + 1) begin
					temp_root_r[i][j] <= temp_root_w[i][j];
				end
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

	function automatic [10:0] shift_poly_1_6;
		input [9:0] i_poly;
		begin
			shift_poly_1_6[0] = i_poly[5];
			shift_poly_1_6[1] = i_poly[0] ^ i_poly[5];	
			shift_poly_1_6[2] = i_poly[1];	
			shift_poly_1_6[3] = i_poly[2];	
			shift_poly_1_6[4] = i_poly[3];	
			shift_poly_1_6[5] = i_poly[4];
			shift_poly_1_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_2_6;
		input [9:0] i_poly;
		begin
			shift_poly_2_6[0] = i_poly[4];
			shift_poly_2_6[1] = i_poly[4] ^ i_poly[5];	
			shift_poly_2_6[2] = i_poly[0] ^ i_poly[5];	
			shift_poly_2_6[3] = i_poly[1];	
			shift_poly_2_6[4] = i_poly[2];	
			shift_poly_2_6[5] = i_poly[3];
			shift_poly_2_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_3_6;
		input [9:0] i_poly;
		begin
			shift_poly_3_6[0] = i_poly[3];
			shift_poly_3_6[1] = i_poly[4] ^ i_poly[3];	
			shift_poly_3_6[2] = i_poly[4] ^ i_poly[5];	
			shift_poly_3_6[3] = i_poly[0] ^ i_poly[5];	
			shift_poly_3_6[4] = i_poly[1];	
			shift_poly_3_6[5] = i_poly[2];
			shift_poly_3_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_4_6;
		input [9:0] i_poly;
		begin
			shift_poly_4_6[0] = i_poly[2];
			shift_poly_4_6[1] = i_poly[2] ^ i_poly[3];	
			shift_poly_4_6[2] = i_poly[4] ^ i_poly[3];	
			shift_poly_4_6[3] = i_poly[4] ^ i_poly[5];	
			shift_poly_4_6[4] = i_poly[0] ^ i_poly[5];	
			shift_poly_4_6[5] = i_poly[1];
			shift_poly_4_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_5_6;
		input [9:0] i_poly;
		begin
			shift_poly_5_6[0] = i_poly[1];
			shift_poly_5_6[1] = i_poly[2] ^ i_poly[1];
			shift_poly_5_6[2] = i_poly[2] ^ i_poly[3];	
			shift_poly_5_6[3] = i_poly[4] ^ i_poly[3];	
			shift_poly_5_6[4] = i_poly[4] ^ i_poly[5];	
			shift_poly_5_6[5] = i_poly[0] ^ i_poly[5];	
			shift_poly_5_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_6_6;
		input [9:0] i_poly;
		begin
			shift_poly_6_6[0] = i_poly[0] ^ i_poly[5];
			shift_poly_6_6[1] = i_poly[0] ^ i_poly[1] ^ i_poly[5];
			shift_poly_6_6[2] = i_poly[2] ^ i_poly[1];
			shift_poly_6_6[3] = i_poly[2] ^ i_poly[3];	
			shift_poly_6_6[4] = i_poly[4] ^ i_poly[3];	
			shift_poly_6_6[5] = i_poly[4] ^ i_poly[5];	
			shift_poly_6_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_7_6;
		input [9:0] i_poly;
		begin
			shift_poly_7_6[0] = i_poly[4] ^ i_poly[5];
			shift_poly_7_6[1] = i_poly[0] ^ i_poly[4];
			shift_poly_7_6[2] = i_poly[0] ^ i_poly[1] ^ i_poly[5];
			shift_poly_7_6[3] = i_poly[2] ^ i_poly[1];	
			shift_poly_7_6[4] = i_poly[2] ^ i_poly[3];	
			shift_poly_7_6[5] = i_poly[4] ^ i_poly[3];	
			shift_poly_7_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_8_6;
		input [9:0] i_poly;
		begin
			shift_poly_8_6[0] = i_poly[4] ^ i_poly[3];
			shift_poly_8_6[1] = i_poly[5] ^ i_poly[3];
			shift_poly_8_6[2] = i_poly[0] ^ i_poly[4];
			shift_poly_8_6[3] = i_poly[1] ^ i_poly[0] ^ i_poly[5];	
			shift_poly_8_6[4] = i_poly[2] ^ i_poly[1];	
			shift_poly_8_6[5] = i_poly[3] ^ i_poly[2];	
			shift_poly_8_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_9_6;
		input [9:0] i_poly;
		begin
			shift_poly_9_6[0] = i_poly[3] ^ i_poly[2];
			shift_poly_9_6[1] = i_poly[4] ^ i_poly[2];
			shift_poly_9_6[2] = i_poly[5] ^ i_poly[3];
			shift_poly_9_6[3] = i_poly[0] ^ i_poly[4];
			shift_poly_9_6[4] = i_poly[1] ^ i_poly[0] ^ i_poly[5];
			shift_poly_9_6[5] = i_poly[2] ^ i_poly[1];
			shift_poly_9_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_10_6;
		input [9:0] i_poly;
		begin
			shift_poly_10_6[0] = i_poly[2] ^ i_poly[1];
			shift_poly_10_6[1] = i_poly[3] ^ i_poly[1];
			shift_poly_10_6[2] = i_poly[4] ^ i_poly[2];
			shift_poly_10_6[3] = i_poly[5] ^ i_poly[3];
			shift_poly_10_6[4] = i_poly[0] ^ i_poly[4];
			shift_poly_10_6[5] = i_poly[1] ^ i_poly[0] ^ i_poly[5];
			shift_poly_10_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_11_6;
		input [9:0] i_poly;
		begin
			shift_poly_11_6[0] = i_poly[1] ^ i_poly[0] ^ i_poly[5];
			shift_poly_11_6[1] = i_poly[2] ^ i_poly[0] ^ i_poly[5];
			shift_poly_11_6[2] = i_poly[3] ^ i_poly[1];
			shift_poly_11_6[3] = i_poly[4] ^ i_poly[2];
			shift_poly_11_6[4] = i_poly[5] ^ i_poly[3];
			shift_poly_11_6[5] = i_poly[0] ^ i_poly[4];
			shift_poly_11_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_12_6;
		input [9:0] i_poly;
		begin
			shift_poly_12_6[0] = i_poly[0] ^ i_poly[4];
			shift_poly_12_6[1] = i_poly[1] ^ i_poly[4] ^ i_poly[5];
			shift_poly_12_6[2] = i_poly[2] ^ i_poly[0] ^ i_poly[5];
			shift_poly_12_6[3] = i_poly[3] ^ i_poly[1];
			shift_poly_12_6[4] = i_poly[4] ^ i_poly[2];
			shift_poly_12_6[5] = i_poly[5] ^ i_poly[3];
			shift_poly_12_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_13_6;
		input [9:0] i_poly;
		begin
			shift_poly_13_6[0] = i_poly[5] ^ i_poly[3];
			shift_poly_13_6[1] = i_poly[0] ^ i_poly[4] ^ i_poly[5] ^ i_poly[3];
			shift_poly_13_6[2] = i_poly[1] ^ i_poly[4] ^ i_poly[5];
			shift_poly_13_6[3] = i_poly[2] ^ i_poly[0] ^ i_poly[5];
			shift_poly_13_6[4] = i_poly[3] ^ i_poly[1];
			shift_poly_13_6[5] = i_poly[4] ^ i_poly[2];
			shift_poly_13_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_14_6;
		input [9:0] i_poly;
		begin
			shift_poly_14_6[0] = i_poly[4] ^ i_poly[2];
			shift_poly_14_6[1] = i_poly[5] ^ i_poly[3] ^ i_poly[4] ^ i_poly[2];
			shift_poly_14_6[2] = i_poly[0] ^ i_poly[4] ^ i_poly[5] ^ i_poly[3];
			shift_poly_14_6[3] = i_poly[1] ^ i_poly[4] ^ i_poly[5];
			shift_poly_14_6[4] = i_poly[2] ^ i_poly[0] ^ i_poly[5];
			shift_poly_14_6[5] = i_poly[3] ^ i_poly[1];
			shift_poly_14_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_15_6;
		input [9:0] i_poly;
		begin
			shift_poly_15_6[0] = i_poly[3] ^ i_poly[1];
			shift_poly_15_6[1] = i_poly[4] ^ i_poly[2] ^ i_poly[3] ^ i_poly[1];
			shift_poly_15_6[2] = i_poly[5] ^ i_poly[3] ^ i_poly[4] ^ i_poly[2];
			shift_poly_15_6[3] = i_poly[0] ^ i_poly[4] ^ i_poly[5] ^ i_poly[3];
			shift_poly_15_6[4] = i_poly[1] ^ i_poly[4] ^ i_poly[5];
			shift_poly_15_6[5] = i_poly[2] ^ i_poly[0] ^ i_poly[5];
			shift_poly_15_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_16_6;
		input [9:0] i_poly;
		begin
			shift_poly_16_6[0] = i_poly[2] ^ i_poly[0] ^ i_poly[5];
			shift_poly_16_6[1] = i_poly[3] ^ i_poly[1] ^ i_poly[2] ^ i_poly[0] ^ i_poly[5];
			shift_poly_16_6[2] = i_poly[4] ^ i_poly[2] ^ i_poly[3] ^ i_poly[1];
			shift_poly_16_6[3] = i_poly[5] ^ i_poly[3] ^ i_poly[4] ^ i_poly[2];
			shift_poly_16_6[4] = i_poly[0] ^ i_poly[4] ^ i_poly[5] ^ i_poly[3];
			shift_poly_16_6[5] = i_poly[1] ^ i_poly[4] ^ i_poly[5];
			shift_poly_16_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_17_6;
		input [9:0] i_poly;
		begin
			shift_poly_17_6[0] = i_poly[1] ^ i_poly[4] ^ i_poly[5];
			shift_poly_17_6[1] = i_poly[2] ^ i_poly[0] ^ i_poly[1] ^ i_poly[4];
			shift_poly_17_6[2] = i_poly[3] ^ i_poly[1] ^ i_poly[2] ^ i_poly[0] ^ i_poly[5];
			shift_poly_17_6[3] = i_poly[4] ^ i_poly[2] ^ i_poly[3] ^ i_poly[1];
			shift_poly_17_6[4] = i_poly[5] ^ i_poly[3] ^ i_poly[4] ^ i_poly[2];
			shift_poly_17_6[5] = i_poly[0] ^ i_poly[4] ^ i_poly[5] ^ i_poly[3];
			shift_poly_17_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_18_6;
		input [9:0] i_poly;
		begin
			shift_poly_18_6[0] = i_poly[0] ^ i_poly[4] ^ i_poly[5] ^ i_poly[3];
			shift_poly_18_6[1] = i_poly[1] ^ i_poly[0] ^ i_poly[3];
			shift_poly_18_6[2] = i_poly[2] ^ i_poly[0] ^ i_poly[1] ^ i_poly[4];
			shift_poly_18_6[3] = i_poly[3] ^ i_poly[1] ^ i_poly[2] ^ i_poly[0] ^ i_poly[5];
			shift_poly_18_6[4] = i_poly[4] ^ i_poly[2] ^ i_poly[3] ^ i_poly[1];
			shift_poly_18_6[5] = i_poly[5] ^ i_poly[3] ^ i_poly[4] ^ i_poly[2];
			shift_poly_18_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_19_6;
		input [9:0] i_poly;
		begin
			shift_poly_19_6[0] = i_poly[5] ^ i_poly[3] ^ i_poly[4] ^ i_poly[2];
			shift_poly_19_6[1] = i_poly[0] ^ i_poly[2];
			shift_poly_19_6[2] = i_poly[1] ^ i_poly[0] ^ i_poly[3];
			shift_poly_19_6[3] = i_poly[2] ^ i_poly[0] ^ i_poly[1] ^ i_poly[4];
			shift_poly_19_6[4] = i_poly[3] ^ i_poly[1] ^ i_poly[2] ^ i_poly[0] ^ i_poly[5];
			shift_poly_19_6[5] = i_poly[4] ^ i_poly[2] ^ i_poly[3] ^ i_poly[1];
			shift_poly_19_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_20_6;
		input [9:0] i_poly;
		begin
			shift_poly_20_6[0] = i_poly[4] ^ i_poly[2] ^ i_poly[3] ^ i_poly[1];
			shift_poly_20_6[1] = i_poly[5] ^ i_poly[1];
			shift_poly_20_6[2] = i_poly[0] ^ i_poly[2];
			shift_poly_20_6[3] = i_poly[1] ^ i_poly[0] ^ i_poly[3];
			shift_poly_20_6[4] = i_poly[2] ^ i_poly[0] ^ i_poly[1] ^ i_poly[4];
			shift_poly_20_6[5] = i_poly[3] ^ i_poly[1] ^ i_poly[2] ^ i_poly[0] ^ i_poly[5];
			shift_poly_20_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_21_6;
		input [9:0] i_poly;
		begin
			shift_poly_21_6[0] = i_poly[3] ^ i_poly[1] ^ i_poly[2] ^ i_poly[0] ^ i_poly[5];
			shift_poly_21_6[1] = i_poly[4] ^ i_poly[0] ^ i_poly[5];
			shift_poly_21_6[2] = i_poly[5] ^ i_poly[1];
			shift_poly_21_6[3] = i_poly[0] ^ i_poly[2];
			shift_poly_21_6[4] = i_poly[1] ^ i_poly[0] ^ i_poly[3];
			shift_poly_21_6[5] = i_poly[2] ^ i_poly[0] ^ i_poly[1] ^ i_poly[4];
			shift_poly_21_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_22_6;
		input [9:0] i_poly;
		begin
			shift_poly_22_6[0] = i_poly[2] ^ i_poly[0] ^ i_poly[1] ^ i_poly[4];
			shift_poly_22_6[1] = i_poly[3] ^ i_poly[5] ^ i_poly[4];
			shift_poly_22_6[2] = i_poly[4] ^ i_poly[0] ^ i_poly[5];
			shift_poly_22_6[3] = i_poly[5] ^ i_poly[1];
			shift_poly_22_6[4] = i_poly[0] ^ i_poly[2];
			shift_poly_22_6[5] = i_poly[1] ^ i_poly[0] ^ i_poly[3];
			shift_poly_22_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_23_6;
		input [9:0] i_poly;
		begin
			shift_poly_23_6[0] = i_poly[1] ^ i_poly[0] ^ i_poly[3];
			shift_poly_23_6[1] = i_poly[2] ^ i_poly[4] ^ i_poly[3];
			shift_poly_23_6[2] = i_poly[3] ^ i_poly[5] ^ i_poly[4];
			shift_poly_23_6[3] = i_poly[4] ^ i_poly[0] ^ i_poly[5];
			shift_poly_23_6[4] = i_poly[5] ^ i_poly[1];
			shift_poly_23_6[5] = i_poly[0] ^ i_poly[2];
			shift_poly_23_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_24_6;
		input [9:0] i_poly;
		begin
			shift_poly_24_6[0] = i_poly[0] ^ i_poly[2];
			shift_poly_24_6[1] = i_poly[1] ^ i_poly[2] ^ i_poly[3];
			shift_poly_24_6[2] = i_poly[2] ^ i_poly[4] ^ i_poly[3];
			shift_poly_24_6[3] = i_poly[3] ^ i_poly[5] ^ i_poly[4];
			shift_poly_24_6[4] = i_poly[4] ^ i_poly[0] ^ i_poly[5];
			shift_poly_24_6[5] = i_poly[5] ^ i_poly[1];
			shift_poly_24_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_25_6;
		input [9:0] i_poly;
		begin
			shift_poly_25_6[0] = i_poly[5] ^ i_poly[1];
			shift_poly_25_6[1] = i_poly[0] ^ i_poly[2] ^ i_poly[5] ^ i_poly[1];
			shift_poly_25_6[2] = i_poly[1] ^ i_poly[2] ^ i_poly[3];
			shift_poly_25_6[3] = i_poly[2] ^ i_poly[4] ^ i_poly[3];
			shift_poly_25_6[4] = i_poly[3] ^ i_poly[5] ^ i_poly[4];
			shift_poly_25_6[5] = i_poly[4] ^ i_poly[0] ^ i_poly[5];
			shift_poly_25_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_26_6;
		input [9:0] i_poly;
		begin
			shift_poly_26_6[0] = i_poly[4] ^ i_poly[0] ^ i_poly[5];
			shift_poly_26_6[1] = i_poly[4] ^ i_poly[1] ^ i_poly[0];
			shift_poly_26_6[2] = i_poly[0] ^ i_poly[2] ^ i_poly[5] ^ i_poly[1];
			shift_poly_26_6[3] = i_poly[1] ^ i_poly[2] ^ i_poly[3];
			shift_poly_26_6[4] = i_poly[2] ^ i_poly[4] ^ i_poly[3];
			shift_poly_26_6[5] = i_poly[3] ^ i_poly[5] ^ i_poly[4];
			shift_poly_26_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_27_6;
		input [9:0] i_poly;
		begin
			shift_poly_27_6[0] = i_poly[3] ^ i_poly[5] ^ i_poly[4];
			shift_poly_27_6[1] = i_poly[3] ^ i_poly[0];
			shift_poly_27_6[2] = i_poly[4] ^ i_poly[1] ^ i_poly[0];
			shift_poly_27_6[3] = i_poly[0] ^ i_poly[2] ^ i_poly[5] ^ i_poly[1];
			shift_poly_27_6[4] = i_poly[1] ^ i_poly[2] ^ i_poly[3];
			shift_poly_27_6[5] = i_poly[2] ^ i_poly[4] ^ i_poly[3];
			shift_poly_27_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_28_6;
		input [9:0] i_poly;
		begin
			shift_poly_28_6[0] = i_poly[2] ^ i_poly[4] ^ i_poly[3];
			shift_poly_28_6[1] = i_poly[2] ^ i_poly[5];
			shift_poly_28_6[2] = i_poly[3] ^ i_poly[0];
			shift_poly_28_6[3] = i_poly[4] ^ i_poly[1] ^ i_poly[0];
			shift_poly_28_6[4] = i_poly[0] ^ i_poly[2] ^ i_poly[5] ^ i_poly[1];
			shift_poly_28_6[5] = i_poly[1] ^ i_poly[2] ^ i_poly[3];
			shift_poly_28_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_29_6;
		input [9:0] i_poly;
		begin
			shift_poly_29_6[0] = i_poly[1] ^ i_poly[2] ^ i_poly[3];
			shift_poly_29_6[1] = i_poly[1] ^ i_poly[4];
			shift_poly_29_6[2] = i_poly[2] ^ i_poly[5];
			shift_poly_29_6[3] = i_poly[3] ^ i_poly[0];
			shift_poly_29_6[4] = i_poly[4] ^ i_poly[1] ^ i_poly[0];
			shift_poly_29_6[5] = i_poly[0] ^ i_poly[2] ^ i_poly[5] ^ i_poly[1];
			shift_poly_29_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_30_6;
		input [9:0] i_poly;
		begin
			shift_poly_30_6[0] = i_poly[0] ^ i_poly[2] ^ i_poly[5] ^ i_poly[1];
			shift_poly_30_6[1] = i_poly[0] ^ i_poly[5] ^ i_poly[3];
			shift_poly_30_6[2] = i_poly[1] ^ i_poly[4];
			shift_poly_30_6[3] = i_poly[2] ^ i_poly[5];
			shift_poly_30_6[4] = i_poly[3] ^ i_poly[0];
			shift_poly_30_6[5] = i_poly[4] ^ i_poly[1] ^ i_poly[0];
			shift_poly_30_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_31_6;
		input [9:0] i_poly;
		begin
			shift_poly_31_6[0] = i_poly[4] ^ i_poly[1] ^ i_poly[0];
			shift_poly_31_6[1] = i_poly[5] ^ i_poly[2] ^ i_poly[4];
			shift_poly_31_6[2] = i_poly[0] ^ i_poly[5] ^ i_poly[3];
			shift_poly_31_6[3] = i_poly[1] ^ i_poly[4];
			shift_poly_31_6[4] = i_poly[2] ^ i_poly[5];
			shift_poly_31_6[5] = i_poly[3] ^ i_poly[0];
			shift_poly_31_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_32_6;
		input [9:0] i_poly;
		begin
			shift_poly_32_6[0] = i_poly[3] ^ i_poly[0];
			shift_poly_32_6[1] = i_poly[4] ^ i_poly[1] ^ i_poly[3];
			shift_poly_32_6[2] = i_poly[5] ^ i_poly[2] ^ i_poly[4];
			shift_poly_32_6[3] = i_poly[0] ^ i_poly[5] ^ i_poly[3];
			shift_poly_32_6[4] = i_poly[1] ^ i_poly[4];
			shift_poly_32_6[5] = i_poly[2] ^ i_poly[5];
			shift_poly_32_6[10:6] = 5'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_1_8;
        input [9:0] i_poly;
        begin
            shift_poly_1_8[0] = i_poly[7];
            shift_poly_1_8[1] = i_poly[0];
            shift_poly_1_8[2] = i_poly[1] ^ i_poly[7];
            shift_poly_1_8[3] = i_poly[2] ^ i_poly[7];
            shift_poly_1_8[4] = i_poly[3] ^ i_poly[7];
            shift_poly_1_8[5] = i_poly[4];
            shift_poly_1_8[6] = i_poly[5];
            shift_poly_1_8[7] = i_poly[6];
            shift_poly_1_8[10:8] = 3'd0;
        end
    endfunction

    function automatic [10:0] shift_poly_2_8;
        input [9:0] i_poly;
        begin
            shift_poly_2_8[0] = i_poly[6];
            shift_poly_2_8[1] = i_poly[7];
            shift_poly_2_8[2] = i_poly[0] ^ i_poly[6];
            shift_poly_2_8[3] = i_poly[1] ^ i_poly[7] ^ i_poly[6];
            shift_poly_2_8[4] = i_poly[2] ^ i_poly[7] ^ i_poly[6];
            shift_poly_2_8[5] = i_poly[3] ^ i_poly[7];
            shift_poly_2_8[6] = i_poly[4];
            shift_poly_2_8[7] = i_poly[5];
            shift_poly_2_8[10:8] = 3'd0;
        end
    endfunction

    function automatic [10:0] shift_poly_3_8;
        input [9:0] i_poly;
        begin
            shift_poly_3_8[0] = i_poly[5];
            shift_poly_3_8[1] = i_poly[6];
            shift_poly_3_8[2] = i_poly[7] ^ i_poly[5];
            shift_poly_3_8[3] = i_poly[0] ^ i_poly[6] ^ i_poly[5];
            shift_poly_3_8[4] = i_poly[1] ^ i_poly[7] ^ i_poly[6] ^ i_poly[5];
            shift_poly_3_8[5] = i_poly[2] ^ i_poly[7] ^ i_poly[6];
            shift_poly_3_8[6] = i_poly[3] ^ i_poly[7];
            shift_poly_3_8[7] = i_poly[4];
            shift_poly_3_8[10:8] = 3'd0;
        end
    endfunction

    function automatic [10:0] shift_poly_4_8;
        input [9:0] i_poly;
        begin
            shift_poly_4_8[0] = i_poly[4];
            shift_poly_4_8[1] = i_poly[5];
            shift_poly_4_8[2] = i_poly[6] ^ i_poly[4];
            shift_poly_4_8[3] = i_poly[7] ^ i_poly[5] ^ i_poly[4];
            shift_poly_4_8[4] = i_poly[0] ^ i_poly[6] ^ i_poly[5] ^ i_poly[4];
            shift_poly_4_8[5] = i_poly[1] ^ i_poly[7] ^ i_poly[6] ^ i_poly[5];
            shift_poly_4_8[6] = i_poly[2] ^ i_poly[7] ^ i_poly[6];
            shift_poly_4_8[7] = i_poly[3] ^ i_poly[7];
            shift_poly_4_8[10:8] = 3'd0;
        end
    endfunction

    function automatic [10:0] shift_poly_5_8;
        input [9:0] i_poly;
        begin
            shift_poly_5_8[0] = i_poly[3] ^ i_poly[7];
            shift_poly_5_8[1] = i_poly[4];
            shift_poly_5_8[2] = i_poly[5] ^ i_poly[3] ^ i_poly[7];
            shift_poly_5_8[3] = i_poly[6] ^ i_poly[4] ^ i_poly[3] ^ i_poly[7];
            shift_poly_5_8[4] = i_poly[5] ^ i_poly[4] ^ i_poly[3];
            shift_poly_5_8[5] = i_poly[0] ^ i_poly[6] ^ i_poly[5] ^ i_poly[4];
            shift_poly_5_8[6] = i_poly[1] ^ i_poly[7] ^ i_poly[6] ^ i_poly[5];
            shift_poly_5_8[7] = i_poly[2] ^ i_poly[7] ^ i_poly[6];
            shift_poly_5_8[10:8] = 3'd0;
        end
    endfunction

    function automatic [10:0] shift_poly_6_8;
        input [9:0] i_poly;
        begin
            shift_poly_6_8[0] = i_poly[2] ^ i_poly[7] ^ i_poly[6];
            shift_poly_6_8[1] = i_poly[3] ^ i_poly[7];
            shift_poly_6_8[2] = i_poly[4] ^ i_poly[2] ^ i_poly[7] ^ i_poly[6];
            shift_poly_6_8[3] = i_poly[5] ^ i_poly[3] ^ i_poly[2] ^ i_poly[6];
            shift_poly_6_8[4] = i_poly[4] ^ i_poly[3] ^ i_poly[2];
            shift_poly_6_8[5] = i_poly[5] ^ i_poly[4] ^ i_poly[3];
            shift_poly_6_8[6] = i_poly[0] ^ i_poly[6] ^ i_poly[5] ^ i_poly[4];
            shift_poly_6_8[7] = i_poly[1] ^ i_poly[7] ^ i_poly[6] ^ i_poly[5];
            shift_poly_6_8[10:8] = 3'd0;
        end
    endfunction

    function automatic [10:0] shift_poly_7_8;
        input [9:0] i_poly;
        begin
            shift_poly_7_8[0] = i_poly[1] ^ i_poly[7] ^ i_poly[6] ^ i_poly[5];
            shift_poly_7_8[1] = i_poly[2] ^ i_poly[7] ^ i_poly[6];
            shift_poly_7_8[2] = i_poly[3] ^ i_poly[1] ^ i_poly[6] ^ i_poly[5];
            shift_poly_7_8[3] = i_poly[4] ^ i_poly[2] ^ i_poly[1] ^ i_poly[5];
            shift_poly_7_8[4] = i_poly[3] ^ i_poly[2] ^ i_poly[1] ^ i_poly[7];
            shift_poly_7_8[5] = i_poly[4] ^ i_poly[3] ^ i_poly[2];
            shift_poly_7_8[6] = i_poly[5] ^ i_poly[4] ^ i_poly[3];
            shift_poly_7_8[7] = i_poly[0] ^ i_poly[6] ^ i_poly[5] ^ i_poly[4];
            shift_poly_7_8[10:8] = 3'd0;
        end
    endfunction

    function automatic [10:0] shift_poly_8_8;
        input [9:0] i_poly;
        begin
            shift_poly_8_8[0] = i_poly[0] ^ i_poly[6] ^ i_poly[5] ^ i_poly[4];
            shift_poly_8_8[1] = i_poly[1] ^ i_poly[7] ^ i_poly[6] ^ i_poly[5];
            shift_poly_8_8[2] = i_poly[2] ^ i_poly[7] ^ i_poly[0] ^ i_poly[5] ^ i_poly[4];
            shift_poly_8_8[3] = i_poly[3] ^ i_poly[1] ^ i_poly[0] ^ i_poly[4];
            shift_poly_8_8[4] = i_poly[2] ^ i_poly[1] ^ i_poly[0] ^ i_poly[6];
            shift_poly_8_8[5] = i_poly[3] ^ i_poly[2] ^ i_poly[1] ^ i_poly[7];
            shift_poly_8_8[6] = i_poly[4] ^ i_poly[3] ^ i_poly[2];
            shift_poly_8_8[7] = i_poly[5] ^ i_poly[4] ^ i_poly[3];
            shift_poly_8_8[10:8] = 3'd0;
        end
    endfunction

	function automatic [10:0] shift_poly_9_8;
        input [9:0] i_poly;
        begin
            shift_poly_9_8[0] = i_poly[5] ^ i_poly[4] ^ i_poly[3];
            shift_poly_9_8[1] = i_poly[0] ^ i_poly[6] ^ i_poly[5] ^ i_poly[4];
            shift_poly_9_8[2] = i_poly[1] ^ i_poly[7] ^ i_poly[6] ^ i_poly[4] ^ i_poly[3];
            shift_poly_9_8[3] = i_poly[2] ^ i_poly[7] ^ i_poly[0] ^ i_poly[3];
            shift_poly_9_8[4] = i_poly[5] ^ i_poly[1] ^ i_poly[0];
            shift_poly_9_8[5] = i_poly[2] ^ i_poly[1] ^ i_poly[0] ^ i_poly[6];
            shift_poly_9_8[6] = i_poly[3] ^ i_poly[2] ^ i_poly[1] ^ i_poly[7];
            shift_poly_9_8[7] = i_poly[4] ^ i_poly[3] ^ i_poly[2];
            shift_poly_9_8[10:8] = 3'd0;
        end
    endfunction

	function automatic [10:0] shift_poly_10_8;
        input [9:0] i_poly;
        begin
            shift_poly_10_8[0] = i_poly[4] ^ i_poly[3] ^ i_poly[2];
            shift_poly_10_8[1] = i_poly[5] ^ i_poly[4] ^ i_poly[3];
            shift_poly_10_8[2] = i_poly[0] ^ i_poly[6] ^ i_poly[5] ^ i_poly[3] ^ i_poly[2];
            shift_poly_10_8[3] = i_poly[1] ^ i_poly[7] ^ i_poly[6] ^ i_poly[2];
            shift_poly_10_8[4] = i_poly[4] ^ i_poly[7] ^ i_poly[0];
            shift_poly_10_8[5] = i_poly[5] ^ i_poly[1] ^ i_poly[0];
            shift_poly_10_8[6] = i_poly[2] ^ i_poly[1] ^ i_poly[0] ^ i_poly[6];
            shift_poly_10_8[7] = i_poly[3] ^ i_poly[2] ^ i_poly[1] ^ i_poly[7];
            shift_poly_10_8[10:8] = 3'd0;
        end
    endfunction

	function automatic [10:0] shift_poly_11_8;
        input [9:0] i_poly;
        begin
            shift_poly_11_8[0] = i_poly[3] ^ i_poly[2] ^ i_poly[1] ^ i_poly[7];
            shift_poly_11_8[1] = i_poly[4] ^ i_poly[3] ^ i_poly[2];
            shift_poly_11_8[2] = i_poly[5] ^ i_poly[4] ^ i_poly[2] ^ i_poly[1] ^ i_poly[7];
            shift_poly_11_8[3] = i_poly[0] ^ i_poly[6] ^ i_poly[5] ^ i_poly[1] ^ i_poly[7];
            shift_poly_11_8[4] = i_poly[3] ^ i_poly[6];
            shift_poly_11_8[5] = i_poly[4] ^ i_poly[7] ^ i_poly[0];
            shift_poly_11_8[6] = i_poly[5] ^ i_poly[1] ^ i_poly[0];
            shift_poly_11_8[7] = i_poly[2] ^ i_poly[1] ^ i_poly[0] ^ i_poly[6];
            shift_poly_11_8[10:8] = 3'd0;
        end
    endfunction

	function automatic [10:0] shift_poly_12_8;
        input [9:0] i_poly;
        begin
            shift_poly_12_8[0] = i_poly[2] ^ i_poly[1] ^ i_poly[0] ^ i_poly[6];
            shift_poly_12_8[1] = i_poly[3] ^ i_poly[2] ^ i_poly[1] ^ i_poly[7];
            shift_poly_12_8[2] = i_poly[4] ^ i_poly[3] ^ i_poly[1] ^ i_poly[0] ^ i_poly[6];
            shift_poly_12_8[3] = i_poly[5] ^ i_poly[4] ^ i_poly[7] ^ i_poly[0] ^ i_poly[6];
            shift_poly_12_8[4] = i_poly[2] ^ i_poly[5] ^ i_poly[7];
            shift_poly_12_8[5] = i_poly[3] ^ i_poly[6];
            shift_poly_12_8[6] = i_poly[4] ^ i_poly[7] ^ i_poly[0];
            shift_poly_12_8[7] = i_poly[5] ^ i_poly[1] ^ i_poly[0];
            shift_poly_12_8[10:8] = 3'd0;
        end
    endfunction

	function automatic [10:0] shift_poly_13_8;
        input [9:0] i_poly;
        begin
            shift_poly_13_8[0] = i_poly[5] ^ i_poly[1] ^ i_poly[0];
            shift_poly_13_8[1] = i_poly[2] ^ i_poly[1] ^ i_poly[0] ^ i_poly[6];
            shift_poly_13_8[2] = i_poly[3] ^ i_poly[2] ^ i_poly[5] ^ i_poly[7] ^ i_poly[0];
            shift_poly_13_8[3] = i_poly[4] ^ i_poly[3] ^ i_poly[5] ^ i_poly[6];
            shift_poly_13_8[4] = i_poly[1] ^ i_poly[4] ^ i_poly[6] ^ i_poly[7];
            shift_poly_13_8[5] = i_poly[2] ^ i_poly[5] ^ i_poly[7];
            shift_poly_13_8[6] = i_poly[3] ^ i_poly[6];
            shift_poly_13_8[7] = i_poly[4] ^ i_poly[7] ^ i_poly[0];
            shift_poly_13_8[10:8] = 3'd0;
        end
    endfunction

	function automatic [10:0] shift_poly_14_8;
        input [9:0] i_poly;
        begin
            shift_poly_14_8[0] = i_poly[4] ^ i_poly[7] ^ i_poly[0];
            shift_poly_14_8[1] = i_poly[5] ^ i_poly[1] ^ i_poly[0];
            shift_poly_14_8[2] = i_poly[2] ^ i_poly[1] ^ i_poly[4] ^ i_poly[6] ^ i_poly[7];
            shift_poly_14_8[3] = i_poly[3] ^ i_poly[2] ^ i_poly[5] ^ i_poly[4];
            shift_poly_14_8[4] = i_poly[0] ^ i_poly[3] ^ i_poly[5] ^ i_poly[6] ^ i_poly[7];
            shift_poly_14_8[5] = i_poly[1] ^ i_poly[4] ^ i_poly[6] ^ i_poly[7];
            shift_poly_14_8[6] = i_poly[2] ^ i_poly[5] ^ i_poly[7];
            shift_poly_14_8[7] = i_poly[3] ^ i_poly[6];
            shift_poly_14_8[10:8] = 3'd0;
        end
    endfunction

	function automatic [10:0] shift_poly_15_8;
        input [9:0] i_poly;
        begin
            shift_poly_15_8[0] = i_poly[3] ^ i_poly[6];
            shift_poly_15_8[1] = i_poly[4] ^ i_poly[7] ^ i_poly[0];
            shift_poly_15_8[2] = i_poly[5] ^ i_poly[1] ^ i_poly[0] ^ i_poly[3] ^ i_poly[6];
            shift_poly_15_8[3] = i_poly[2] ^ i_poly[1] ^ i_poly[4] ^ i_poly[3] ^ i_poly[7];
            shift_poly_15_8[4] = i_poly[6] ^ i_poly[2] ^ i_poly[5] ^ i_poly[4];
            shift_poly_15_8[5] = i_poly[0] ^ i_poly[3] ^ i_poly[5] ^ i_poly[6] ^ i_poly[7];
            shift_poly_15_8[6] = i_poly[1] ^ i_poly[4] ^ i_poly[6] ^ i_poly[7];
            shift_poly_15_8[7] = i_poly[2] ^ i_poly[5] ^ i_poly[7];
            shift_poly_15_8[10:8] = 3'd0;
        end
    endfunction

	function automatic [10:0] shift_poly_16_8;
        input [9:0] i_poly;
        begin
            shift_poly_16_8[0] = i_poly[2] ^ i_poly[5] ^ i_poly[7];
            shift_poly_16_8[1] = i_poly[3] ^ i_poly[6];
            shift_poly_16_8[2] = i_poly[4] ^ i_poly[2] ^ i_poly[0] ^ i_poly[5];
            shift_poly_16_8[3] = i_poly[2] ^ i_poly[1] ^ i_poly[0] ^ i_poly[3] ^ i_poly[6] ^ i_poly[7];
            shift_poly_16_8[4] = i_poly[5] ^ i_poly[1] ^ i_poly[4] ^ i_poly[3];
            shift_poly_16_8[5] = i_poly[6] ^ i_poly[2] ^ i_poly[5] ^ i_poly[4];
            shift_poly_16_8[6] = i_poly[0] ^ i_poly[3] ^ i_poly[5] ^ i_poly[6] ^ i_poly[7];
            shift_poly_16_8[7] = i_poly[1] ^ i_poly[4] ^ i_poly[6] ^ i_poly[7];
            shift_poly_16_8[10:8] = 3'd0;
        end
    endfunction

	function automatic [10:0] shift_poly_17_8;
        input [9:0] i_poly;
        begin
            shift_poly_17_8[0] = i_poly[1] ^ i_poly[4] ^ i_poly[6] ^ i_poly[7];
            shift_poly_17_8[1] = i_poly[2] ^ i_poly[5] ^ i_poly[7];
            shift_poly_17_8[2] = i_poly[3] ^ i_poly[7] ^ i_poly[1] ^ i_poly[4];
            shift_poly_17_8[3] = i_poly[6] ^ i_poly[2] ^ i_poly[0] ^ i_poly[5] ^ i_poly[1] ^ i_poly[7];
            shift_poly_17_8[4] = i_poly[2] ^ i_poly[4] ^ i_poly[0] ^ i_poly[3];
            shift_poly_17_8[5] = i_poly[5] ^ i_poly[1] ^ i_poly[4] ^ i_poly[3];
            shift_poly_17_8[6] = i_poly[6] ^ i_poly[2] ^ i_poly[5] ^ i_poly[4];
            shift_poly_17_8[7] = i_poly[0] ^ i_poly[3] ^ i_poly[5] ^ i_poly[6] ^ i_poly[7];
            shift_poly_17_8[10:8] = 3'd0;
        end
    endfunction

	function automatic [10:0] shift_poly_18_8;
        input [9:0] i_poly;
        begin
            shift_poly_18_8[0] = i_poly[0] ^ i_poly[3] ^ i_poly[5] ^ i_poly[6] ^ i_poly[7];
            shift_poly_18_8[1] = i_poly[1] ^ i_poly[4] ^ i_poly[6] ^ i_poly[7];
            shift_poly_18_8[2] = i_poly[2] ^ i_poly[3] ^ i_poly[6] ^ i_poly[0];
            shift_poly_18_8[3] = i_poly[0] ^ i_poly[5] ^ i_poly[1] ^ i_poly[4] ^ i_poly[6];
            shift_poly_18_8[4] = i_poly[3] ^ i_poly[2] ^ i_poly[1];
            shift_poly_18_8[5] = i_poly[2] ^ i_poly[4] ^ i_poly[0] ^ i_poly[3];
            shift_poly_18_8[6] = i_poly[5] ^ i_poly[1] ^ i_poly[4] ^ i_poly[3];
            shift_poly_18_8[7] = i_poly[6] ^ i_poly[2] ^ i_poly[5] ^ i_poly[4];
            shift_poly_18_8[10:8] = 3'd0;
        end
    endfunction

	function automatic [10:0] shift_poly_19_8;
        input [9:0] i_poly;
        begin
            shift_poly_19_8[0] = i_poly[6] ^ i_poly[2] ^ i_poly[5] ^ i_poly[4];
            shift_poly_19_8[1] = i_poly[0] ^ i_poly[3] ^ i_poly[5] ^ i_poly[6] ^ i_poly[7];
            shift_poly_19_8[2] = i_poly[1] ^ i_poly[2] ^ i_poly[5] ^ i_poly[7];
            shift_poly_19_8[3] = i_poly[5] ^ i_poly[3] ^ i_poly[4] ^ i_poly[0];
            shift_poly_19_8[4] = i_poly[0] ^ i_poly[2] ^ i_poly[1];
            shift_poly_19_8[5] = i_poly[3] ^ i_poly[2] ^ i_poly[1];
            shift_poly_19_8[6] = i_poly[2] ^ i_poly[4] ^ i_poly[0] ^ i_poly[3];
            shift_poly_19_8[7] = i_poly[5] ^ i_poly[1] ^ i_poly[4] ^ i_poly[3];
            shift_poly_19_8[10:8] = 3'd0;
        end
    endfunction

	function automatic [10:0] shift_poly_20_8;
        input [9:0] i_poly;
        begin
            shift_poly_20_8[0] = i_poly[5] ^ i_poly[1] ^ i_poly[4] ^ i_poly[3];
            shift_poly_20_8[1] = i_poly[6] ^ i_poly[2] ^ i_poly[5] ^ i_poly[4];
            shift_poly_20_8[2] = i_poly[0] ^ i_poly[1] ^ i_poly[4] ^ i_poly[6] ^ i_poly[7];
            shift_poly_20_8[3] = i_poly[4] ^ i_poly[2] ^ i_poly[3] ^ i_poly[7];
            shift_poly_20_8[4] = i_poly[0] ^ i_poly[1];
            shift_poly_20_8[5] = i_poly[0] ^ i_poly[2] ^ i_poly[1];
            shift_poly_20_8[6] = i_poly[3] ^ i_poly[2] ^ i_poly[1];
            shift_poly_20_8[7] = i_poly[2] ^ i_poly[4] ^ i_poly[0] ^ i_poly[3];
            shift_poly_20_8[10:8] = 3'd0;
        end
    endfunction

	function automatic [10:0] shift_poly_21_8;
        input [9:0] i_poly;
        begin
            shift_poly_21_8[0] = i_poly[2] ^ i_poly[4] ^ i_poly[0] ^ i_poly[3];
            shift_poly_21_8[1] = i_poly[5] ^ i_poly[1] ^ i_poly[4] ^ i_poly[3];
            shift_poly_21_8[2] = i_poly[6] ^ i_poly[0] ^ i_poly[5] ^ i_poly[3];
            shift_poly_21_8[3] = i_poly[2] ^ i_poly[1] ^ i_poly[3] ^ i_poly[6] ^ i_poly[7];
            shift_poly_21_8[4] = i_poly[0] ^ i_poly[7];
            shift_poly_21_8[5] = i_poly[0] ^ i_poly[1];
            shift_poly_21_8[6] = i_poly[0] ^ i_poly[2] ^ i_poly[1];
            shift_poly_21_8[7] = i_poly[3] ^ i_poly[2] ^ i_poly[1];
            shift_poly_21_8[10:8] = 3'd0;
        end
    endfunction

	function automatic [10:0] shift_poly_22_8;
        input [9:0] i_poly;
        begin
            shift_poly_22_8[0] = i_poly[3] ^ i_poly[2] ^ i_poly[1];
            shift_poly_22_8[1] = i_poly[2] ^ i_poly[4] ^ i_poly[0] ^ i_poly[3];
            shift_poly_22_8[2] = i_poly[5] ^ i_poly[2] ^ i_poly[4];
            shift_poly_22_8[3] = i_poly[6] ^ i_poly[0] ^ i_poly[5] ^ i_poly[2] ^ i_poly[1];
            shift_poly_22_8[4] = i_poly[6] ^ i_poly[7];
            shift_poly_22_8[5] = i_poly[0] ^ i_poly[7];
            shift_poly_22_8[6] = i_poly[0] ^ i_poly[1];
            shift_poly_22_8[7] = i_poly[0] ^ i_poly[2] ^ i_poly[1];
            shift_poly_22_8[10:8] = 3'd0;
        end
    endfunction

	function automatic [10:0] shift_poly_23_8;
        input [9:0] i_poly;
        begin
            shift_poly_23_8[0] = i_poly[0] ^ i_poly[2] ^ i_poly[1];
            shift_poly_23_8[1] = i_poly[3] ^ i_poly[2] ^ i_poly[1];
            shift_poly_23_8[2] = i_poly[1] ^ i_poly[4] ^ i_poly[3];
            shift_poly_23_8[3] = i_poly[5] ^ i_poly[0] ^ i_poly[4] ^ i_poly[1];
            shift_poly_23_8[4] = i_poly[6] ^ i_poly[5];
            shift_poly_23_8[5] = i_poly[6] ^ i_poly[7];
            shift_poly_23_8[6] = i_poly[0] ^ i_poly[7];
            shift_poly_23_8[7] = i_poly[0] ^ i_poly[1];
            shift_poly_23_8[10:8] = 3'd0;
        end
    endfunction

	function automatic [10:0] shift_poly_24_8;
        input [9:0] i_poly;
        begin
            shift_poly_24_8[0] = i_poly[0] ^ i_poly[1];
            shift_poly_24_8[1] = i_poly[0] ^ i_poly[2] ^ i_poly[1];
            shift_poly_24_8[2] = i_poly[3] ^ i_poly[2] ^ i_poly[0];
            shift_poly_24_8[3] = i_poly[0] ^ i_poly[4] ^ i_poly[3];
            shift_poly_24_8[4] = i_poly[5] ^ i_poly[4];
            shift_poly_24_8[5] = i_poly[6] ^ i_poly[5];
            shift_poly_24_8[6] = i_poly[6] ^ i_poly[7];
            shift_poly_24_8[7] = i_poly[0] ^ i_poly[7];
            shift_poly_24_8[10:8] = 3'd0;
        end
    endfunction

	function automatic [10:0] shift_poly_25_8;
		input [9:0] i_poly;
		begin
			shift_poly_25_8[0] = i_poly[0] ^ i_poly[7];
			shift_poly_25_8[1] = i_poly[0] ^ i_poly[1];
			shift_poly_25_8[2] = i_poly[7] ^ i_poly[2] ^ i_poly[1];
			shift_poly_25_8[3] = i_poly[3] ^ i_poly[2] ^ i_poly[7];
			shift_poly_25_8[4] = i_poly[7] ^ i_poly[4] ^ i_poly[3];
			shift_poly_25_8[5] = i_poly[5] ^ i_poly[4];
			shift_poly_25_8[6] = i_poly[6] ^ i_poly[5];
			shift_poly_25_8[7] = i_poly[6] ^ i_poly[7];
			shift_poly_25_8[10:8] = 3'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_26_8;
		input [9:0] i_poly;
		begin
			shift_poly_26_8[0] = i_poly[6] ^ i_poly[7];
			shift_poly_26_8[1] = i_poly[0] ^ i_poly[7];
			shift_poly_26_8[2] = i_poly[0] ^ i_poly[1] ^ i_poly[6] ^ i_poly[7];
			shift_poly_26_8[3] = i_poly[6] ^ i_poly[2] ^ i_poly[1];
			shift_poly_26_8[4] = i_poly[3] ^ i_poly[2] ^ i_poly[6];
			shift_poly_26_8[5] = i_poly[7] ^ i_poly[4] ^ i_poly[3];
			shift_poly_26_8[6] = i_poly[5] ^ i_poly[4];
			shift_poly_26_8[7] = i_poly[6] ^ i_poly[5];
			shift_poly_26_8[10:8] = 3'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_27_8;
		input [9:0] i_poly;
		begin
			shift_poly_27_8[0] = i_poly[6] ^ i_poly[5];
			shift_poly_27_8[1] = i_poly[6] ^ i_poly[7];
			shift_poly_27_8[2] = i_poly[0] ^ i_poly[7] ^ i_poly[6] ^ i_poly[5];
			shift_poly_27_8[3] = i_poly[0] ^ i_poly[1] ^ i_poly[5] ^ i_poly[7];
			shift_poly_27_8[4] = i_poly[5] ^ i_poly[2] ^ i_poly[1];
			shift_poly_27_8[5] = i_poly[3] ^ i_poly[2] ^ i_poly[6];
			shift_poly_27_8[6] = i_poly[7] ^ i_poly[4] ^ i_poly[3];
			shift_poly_27_8[7] = i_poly[5] ^ i_poly[4];
			shift_poly_27_8[10:8] = 3'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_28_8;
		input [9:0] i_poly;
		begin
			shift_poly_28_8[0] = i_poly[5] ^ i_poly[4];
			shift_poly_28_8[1] = i_poly[6] ^ i_poly[5];
			shift_poly_28_8[2] = i_poly[6] ^ i_poly[7] ^ i_poly[5] ^ i_poly[4];
			shift_poly_28_8[3] = i_poly[0] ^ i_poly[7] ^ i_poly[6] ^ i_poly[4];
			shift_poly_28_8[4] = i_poly[0] ^ i_poly[1] ^ i_poly[4] ^ i_poly[7];
			shift_poly_28_8[5] = i_poly[5] ^ i_poly[2] ^ i_poly[1];
			shift_poly_28_8[6] = i_poly[3] ^ i_poly[2] ^ i_poly[6];
			shift_poly_28_8[7] = i_poly[7] ^ i_poly[4] ^ i_poly[3];
			shift_poly_28_8[10:8] = 3'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_29_8;
		input [9:0] i_poly;
		begin
			shift_poly_29_8[0] = i_poly[7] ^ i_poly[4] ^ i_poly[3];
			shift_poly_29_8[1] = i_poly[5] ^ i_poly[4];
			shift_poly_29_8[2] = i_poly[6] ^ i_poly[5] ^ i_poly[7] ^ i_poly[4] ^ i_poly[3];
			shift_poly_29_8[3] = i_poly[6] ^ i_poly[3] ^ i_poly[5];
			shift_poly_29_8[4] = i_poly[0] ^ i_poly[3] ^ i_poly[6];
			shift_poly_29_8[5] = i_poly[0] ^ i_poly[1] ^ i_poly[4] ^ i_poly[7];
			shift_poly_29_8[6] = i_poly[5] ^ i_poly[2] ^ i_poly[1];
			shift_poly_29_8[7] = i_poly[3] ^ i_poly[2] ^ i_poly[6];
			shift_poly_29_8[10:8] = 3'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_30_8;
		input [9:0] i_poly;
		begin
			shift_poly_30_8[0] = i_poly[3] ^ i_poly[2] ^ i_poly[6];
			shift_poly_30_8[1] = i_poly[7] ^ i_poly[4] ^ i_poly[3];
			shift_poly_30_8[2] = i_poly[5] ^ i_poly[4] ^ i_poly[3] ^ i_poly[2] ^ i_poly[6];
			shift_poly_30_8[3] = i_poly[2] ^ i_poly[5] ^ i_poly[7] ^ i_poly[4];
			shift_poly_30_8[4] = i_poly[2] ^ i_poly[5];
			shift_poly_30_8[5] = i_poly[0] ^ i_poly[3] ^ i_poly[6];
			shift_poly_30_8[6] = i_poly[0] ^ i_poly[1] ^ i_poly[4] ^ i_poly[7];
			shift_poly_30_8[7] = i_poly[5] ^ i_poly[2] ^ i_poly[1];
			shift_poly_30_8[10:8] = 3'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_31_8;
		input [9:0] i_poly;
		begin
			shift_poly_31_8[0] = i_poly[5] ^ i_poly[2] ^ i_poly[1];
			shift_poly_31_8[1] = i_poly[3] ^ i_poly[2] ^ i_poly[6];
			shift_poly_31_8[2] = i_poly[7] ^ i_poly[4] ^ i_poly[3] ^ i_poly[5] ^ i_poly[2] ^ i_poly[1];
			shift_poly_31_8[3] = i_poly[1] ^ i_poly[4] ^ i_poly[3] ^ i_poly[6];
			shift_poly_31_8[4] = i_poly[1] ^ i_poly[4] ^ i_poly[7];
			shift_poly_31_8[5] = i_poly[2] ^ i_poly[5];
			shift_poly_31_8[6] = i_poly[0] ^ i_poly[3] ^ i_poly[6];
			shift_poly_31_8[7] = i_poly[0] ^ i_poly[1] ^ i_poly[4] ^ i_poly[7];
			shift_poly_31_8[10:8] = 3'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_32_8;
		input [9:0] i_poly;
		begin
			shift_poly_32_8[0] = i_poly[0] ^ i_poly[1] ^ i_poly[4] ^ i_poly[7];
			shift_poly_32_8[1] = i_poly[5] ^ i_poly[2] ^ i_poly[1];
			shift_poly_32_8[2] = i_poly[3] ^ i_poly[2] ^ i_poly[6] ^ i_poly[0] ^ i_poly[1] ^ i_poly[4] ^ i_poly[7];
			shift_poly_32_8[3] = i_poly[0] ^ i_poly[2] ^ i_poly[3] ^ i_poly[5];
			shift_poly_32_8[4] = i_poly[0] ^ i_poly[7] ^ i_poly[3] ^ i_poly[6];
			shift_poly_32_8[5] = i_poly[1] ^ i_poly[4] ^ i_poly[7];
			shift_poly_32_8[6] = i_poly[2] ^ i_poly[5];
			shift_poly_32_8[7] = i_poly[0] ^ i_poly[3] ^ i_poly[6];
			shift_poly_32_8[10:8] = 3'd0;
		end
	endfunction

	function automatic [10:0] shift_poly_1_10;
        input [9:0] i_poly;
        begin
            shift_poly_1_10[0] = i_poly[9];
            shift_poly_1_10[1] = i_poly[0];
            shift_poly_1_10[2] = i_poly[1];
            shift_poly_1_10[3] = i_poly[2] ^ i_poly[9];
            shift_poly_1_10[4] = i_poly[3];
            shift_poly_1_10[5] = i_poly[4];
            shift_poly_1_10[6] = i_poly[5];
            shift_poly_1_10[7] = i_poly[6];
            shift_poly_1_10[8] = i_poly[7];
            shift_poly_1_10[9] = i_poly[8];
            shift_poly_1_10[10] = 1'b0;
        end
    endfunction

    function automatic [10:0] shift_poly_2_10;
        input [9:0] i_poly;
        begin
            shift_poly_2_10[0] = i_poly[8];
            shift_poly_2_10[1] = i_poly[9];
            shift_poly_2_10[2] = i_poly[0];
            shift_poly_2_10[3] = i_poly[1] ^ i_poly[8];
            shift_poly_2_10[4] = i_poly[2] ^ i_poly[9];
            shift_poly_2_10[5] = i_poly[3];
            shift_poly_2_10[6] = i_poly[4];
            shift_poly_2_10[7] = i_poly[5];
            shift_poly_2_10[8] = i_poly[6];
            shift_poly_2_10[9] = i_poly[7];
            shift_poly_2_10[10] = 1'b0;
        end
    endfunction

    function automatic [10:0] shift_poly_3_10;
        input [9:0] i_poly;
        begin
            shift_poly_3_10[0] = i_poly[7];
            shift_poly_3_10[1] = i_poly[8];
            shift_poly_3_10[2] = i_poly[9];
            shift_poly_3_10[3] = i_poly[0] ^ i_poly[7];
            shift_poly_3_10[4] = i_poly[1] ^ i_poly[8];
            shift_poly_3_10[5] = i_poly[2] ^ i_poly[9];
            shift_poly_3_10[6] = i_poly[3];
            shift_poly_3_10[7] = i_poly[4];
            shift_poly_3_10[8] = i_poly[5];
            shift_poly_3_10[9] = i_poly[6];
            shift_poly_3_10[10] = 1'b0;
        end
    endfunction

    function automatic [10:0] shift_poly_4_10;
        input [9:0] i_poly;
        begin
            shift_poly_4_10[0] = i_poly[6];
            shift_poly_4_10[1] = i_poly[7];
            shift_poly_4_10[2] = i_poly[8];
            shift_poly_4_10[3] = i_poly[9] ^ i_poly[6];
            shift_poly_4_10[4] = i_poly[0] ^ i_poly[7];
            shift_poly_4_10[5] = i_poly[1] ^ i_poly[8];
            shift_poly_4_10[6] = i_poly[2] ^ i_poly[9];
            shift_poly_4_10[7] = i_poly[3];
            shift_poly_4_10[8] = i_poly[4];
            shift_poly_4_10[9] = i_poly[5];
            shift_poly_4_10[10] = 1'b0;
        end
    endfunction

    function automatic [10:0] shift_poly_5_10;
        input [9:0] i_poly;
        begin
            shift_poly_5_10[0] = i_poly[5];
            shift_poly_5_10[1] = i_poly[6];
            shift_poly_5_10[2] = i_poly[7];
            shift_poly_5_10[3] = i_poly[8] ^ i_poly[5];
            shift_poly_5_10[4] = i_poly[9] ^ i_poly[6];
            shift_poly_5_10[5] = i_poly[0] ^ i_poly[7];
            shift_poly_5_10[6] = i_poly[1] ^ i_poly[8];
            shift_poly_5_10[7] = i_poly[2] ^ i_poly[9];
            shift_poly_5_10[8] = i_poly[3];
            shift_poly_5_10[9] = i_poly[4];
            shift_poly_5_10[10] = 1'b0;
        end
    endfunction

    function automatic [10:0] shift_poly_6_10;
        input [9:0] i_poly;
        begin
            shift_poly_6_10[0] = i_poly[4];
            shift_poly_6_10[1] = i_poly[5];
            shift_poly_6_10[2] = i_poly[6];
            shift_poly_6_10[3] = i_poly[7] ^ i_poly[4];
            shift_poly_6_10[4] = i_poly[8] ^ i_poly[5];
            shift_poly_6_10[5] = i_poly[9] ^ i_poly[6];
            shift_poly_6_10[6] = i_poly[0] ^ i_poly[7];
            shift_poly_6_10[7] = i_poly[1] ^ i_poly[8];
            shift_poly_6_10[8] = i_poly[2] ^ i_poly[9];
            shift_poly_6_10[9] = i_poly[3];
            shift_poly_6_10[10] = 1'b0;
        end
    endfunction

    function automatic [10:0] shift_poly_7_10;
        input [9:0] i_poly;
        begin
            shift_poly_7_10[0] = i_poly[3];
            shift_poly_7_10[1] = i_poly[4];
            shift_poly_7_10[2] = i_poly[5];
            shift_poly_7_10[3] = i_poly[6] ^ i_poly[3];
            shift_poly_7_10[4] = i_poly[7] ^ i_poly[4];
            shift_poly_7_10[5] = i_poly[8] ^ i_poly[5];
            shift_poly_7_10[6] = i_poly[9] ^ i_poly[6];
            shift_poly_7_10[7] = i_poly[0] ^ i_poly[7];
            shift_poly_7_10[8] = i_poly[1] ^ i_poly[8];
            shift_poly_7_10[9] = i_poly[2] ^ i_poly[9];
            shift_poly_7_10[10] = 1'b0;
        end
    endfunction

    function automatic [10:0] shift_poly_8_10;
        input [9:0] i_poly;
        begin
            shift_poly_8_10[0] = i_poly[2] ^ i_poly[9];
            shift_poly_8_10[1] = i_poly[3];
            shift_poly_8_10[2] = i_poly[4];
            shift_poly_8_10[3] = i_poly[5] ^ i_poly[9] ^ i_poly[2];
            shift_poly_8_10[4] = i_poly[6] ^ i_poly[3];
            shift_poly_8_10[5] = i_poly[7] ^ i_poly[4];
            shift_poly_8_10[6] = i_poly[8] ^ i_poly[5];
            shift_poly_8_10[7] = i_poly[9] ^ i_poly[6];
            shift_poly_8_10[8] = i_poly[0] ^ i_poly[7];
            shift_poly_8_10[9] = i_poly[1] ^ i_poly[8];
            shift_poly_8_10[10] = 1'b0;
        end
    endfunction

	function automatic [10:0] shift_poly_9_10;
        input [9:0] i_poly;
        begin
            shift_poly_9_10[0] = i_poly[1] ^ i_poly[8];
            shift_poly_9_10[1] = i_poly[2] ^ i_poly[9];
            shift_poly_9_10[2] = i_poly[3];
            shift_poly_9_10[3] = i_poly[4] ^ i_poly[1] ^ i_poly[8];
            shift_poly_9_10[4] = i_poly[5] ^ i_poly[9] ^ i_poly[2];
            shift_poly_9_10[5] = i_poly[6] ^ i_poly[3];
            shift_poly_9_10[6] = i_poly[7] ^ i_poly[4];
            shift_poly_9_10[7] = i_poly[8] ^ i_poly[5];
            shift_poly_9_10[8] = i_poly[9] ^ i_poly[6];
            shift_poly_9_10[9] = i_poly[0] ^ i_poly[7];
            shift_poly_9_10[10] = 1'b0;
        end
    endfunction

	function automatic [10:0] shift_poly_10_10;
        input [9:0] i_poly;
        begin
            shift_poly_10_10[0] = i_poly[0] ^ i_poly[7];
            shift_poly_10_10[1] = i_poly[1] ^ i_poly[8];
            shift_poly_10_10[2] = i_poly[2] ^ i_poly[9];
            shift_poly_10_10[3] = i_poly[3] ^ i_poly[0] ^ i_poly[7];
            shift_poly_10_10[4] = i_poly[4] ^ i_poly[1] ^ i_poly[8];
            shift_poly_10_10[5] = i_poly[5] ^ i_poly[9] ^ i_poly[2];
            shift_poly_10_10[6] = i_poly[6] ^ i_poly[3];
            shift_poly_10_10[7] = i_poly[7] ^ i_poly[4];
            shift_poly_10_10[8] = i_poly[8] ^ i_poly[5];
            shift_poly_10_10[9] = i_poly[9] ^ i_poly[6];
            shift_poly_10_10[10] = 1'b0;
        end
    endfunction

	function automatic [10:0] shift_poly_11_10;
        input [9:0] i_poly;
        begin
            shift_poly_11_10[0] = i_poly[9] ^ i_poly[6];
            shift_poly_11_10[1] = i_poly[0] ^ i_poly[7];
            shift_poly_11_10[2] = i_poly[1] ^ i_poly[8];
            shift_poly_11_10[3] = i_poly[2] ^ i_poly[6];
            shift_poly_11_10[4] = i_poly[3] ^ i_poly[0] ^ i_poly[7];
            shift_poly_11_10[5] = i_poly[4] ^ i_poly[1] ^ i_poly[8];
            shift_poly_11_10[6] = i_poly[5] ^ i_poly[9] ^ i_poly[2];
            shift_poly_11_10[7] = i_poly[6] ^ i_poly[3];
            shift_poly_11_10[8] = i_poly[7] ^ i_poly[4];
            shift_poly_11_10[9] = i_poly[8] ^ i_poly[5];
            shift_poly_11_10[10] = 1'b0;
        end
    endfunction

	function automatic [10:0] shift_poly_12_10;
        input [9:0] i_poly;
        begin
            shift_poly_12_10[0] = i_poly[8] ^ i_poly[5];
            shift_poly_12_10[1] = i_poly[9] ^ i_poly[6];
            shift_poly_12_10[2] = i_poly[0] ^ i_poly[7];
            shift_poly_12_10[3] = i_poly[1] ^ i_poly[5];
            shift_poly_12_10[4] = i_poly[2] ^ i_poly[6];
            shift_poly_12_10[5] = i_poly[3] ^ i_poly[0] ^ i_poly[7];
            shift_poly_12_10[6] = i_poly[4] ^ i_poly[1] ^ i_poly[8];
            shift_poly_12_10[7] = i_poly[5] ^ i_poly[9] ^ i_poly[2];
            shift_poly_12_10[8] = i_poly[6] ^ i_poly[3];
            shift_poly_12_10[9] = i_poly[7] ^ i_poly[4];
            shift_poly_12_10[10] = 1'b0;
        end
    endfunction

	function automatic [10:0] shift_poly_13_10;
		input [9:0] i_poly;
		begin
			shift_poly_13_10[0] = i_poly[7] ^ i_poly[4];
			shift_poly_13_10[1] = i_poly[8] ^ i_poly[5];
			shift_poly_13_10[2] = i_poly[9] ^ i_poly[6];
			shift_poly_13_10[3] = i_poly[0] ^ i_poly[4];
			shift_poly_13_10[4] = i_poly[1] ^ i_poly[5];
			shift_poly_13_10[5] = i_poly[2] ^ i_poly[6];
			shift_poly_13_10[6] = i_poly[3] ^ i_poly[0] ^ i_poly[7];
			shift_poly_13_10[7] = i_poly[4] ^ i_poly[1] ^ i_poly[8];
			shift_poly_13_10[8] = i_poly[5] ^ i_poly[9] ^ i_poly[2];
			shift_poly_13_10[9] = i_poly[6] ^ i_poly[3];
			shift_poly_13_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_14_10;
		input [9:0] i_poly;
		begin
			shift_poly_14_10[0] = i_poly[6] ^ i_poly[3];
			shift_poly_14_10[1] = i_poly[7] ^ i_poly[4];
			shift_poly_14_10[2] = i_poly[8] ^ i_poly[5];
			shift_poly_14_10[3] = i_poly[9] ^ i_poly[3];
			shift_poly_14_10[4] = i_poly[0] ^ i_poly[4];
			shift_poly_14_10[5] = i_poly[1] ^ i_poly[5];
			shift_poly_14_10[6] = i_poly[2] ^ i_poly[6];
			shift_poly_14_10[7] = i_poly[3] ^ i_poly[0] ^ i_poly[7];
			shift_poly_14_10[8] = i_poly[4] ^ i_poly[1] ^ i_poly[8];
			shift_poly_14_10[9] = i_poly[5] ^ i_poly[9] ^ i_poly[2];
			shift_poly_14_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_15_10;
		input [9:0] i_poly;
		begin
			shift_poly_15_10[0] = i_poly[5] ^ i_poly[9] ^ i_poly[2];
			shift_poly_15_10[1] = i_poly[6] ^ i_poly[3];
			shift_poly_15_10[2] = i_poly[7] ^ i_poly[4];
			shift_poly_15_10[3] = i_poly[8] ^ i_poly[2] ^ i_poly[9];
			shift_poly_15_10[4] = i_poly[9] ^ i_poly[3];
			shift_poly_15_10[5] = i_poly[0] ^ i_poly[4];
			shift_poly_15_10[6] = i_poly[1] ^ i_poly[5];
			shift_poly_15_10[7] = i_poly[2] ^ i_poly[6];
			shift_poly_15_10[8] = i_poly[3] ^ i_poly[0] ^ i_poly[7];
			shift_poly_15_10[9] = i_poly[4] ^ i_poly[1] ^ i_poly[8];
			shift_poly_15_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_16_10;
		input [9:0] i_poly;
		begin
			shift_poly_16_10[0] = i_poly[4] ^ i_poly[1] ^ i_poly[8];
			shift_poly_16_10[1] = i_poly[5] ^ i_poly[9] ^ i_poly[2];
			shift_poly_16_10[2] = i_poly[6] ^ i_poly[3];
			shift_poly_16_10[3] = i_poly[7] ^ i_poly[1] ^ i_poly[8];
			shift_poly_16_10[4] = i_poly[8] ^ i_poly[2] ^ i_poly[9];
			shift_poly_16_10[5] = i_poly[9] ^ i_poly[3];
			shift_poly_16_10[6] = i_poly[0] ^ i_poly[4];
			shift_poly_16_10[7] = i_poly[1] ^ i_poly[5];
			shift_poly_16_10[8] = i_poly[2] ^ i_poly[6];
			shift_poly_16_10[9] = i_poly[3] ^ i_poly[0] ^ i_poly[7];
			shift_poly_16_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_17_10;
		input [9:0] i_poly;
		begin
			shift_poly_17_10[0] = i_poly[3] ^ i_poly[0] ^ i_poly[7];
			shift_poly_17_10[1] = i_poly[4] ^ i_poly[1] ^ i_poly[8];
			shift_poly_17_10[2] = i_poly[5] ^ i_poly[9] ^ i_poly[2];
			shift_poly_17_10[3] = i_poly[6] ^ i_poly[0] ^ i_poly[7];
			shift_poly_17_10[4] = i_poly[7] ^ i_poly[1] ^ i_poly[8];
			shift_poly_17_10[5] = i_poly[8] ^ i_poly[2] ^ i_poly[9];
			shift_poly_17_10[6] = i_poly[9] ^ i_poly[3];
			shift_poly_17_10[7] = i_poly[0] ^ i_poly[4];
			shift_poly_17_10[8] = i_poly[1] ^ i_poly[5];
			shift_poly_17_10[9] = i_poly[2] ^ i_poly[6];
			shift_poly_17_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_18_10;
		input [9:0] i_poly;
		begin
			shift_poly_18_10[0] = i_poly[2] ^ i_poly[6];
			shift_poly_18_10[1] = i_poly[3] ^ i_poly[0] ^ i_poly[7];
			shift_poly_18_10[2] = i_poly[4] ^ i_poly[1] ^ i_poly[8];
			shift_poly_18_10[3] = i_poly[5] ^ i_poly[9] ^ i_poly[6];
			shift_poly_18_10[4] = i_poly[6] ^ i_poly[0] ^ i_poly[7];
			shift_poly_18_10[5] = i_poly[7] ^ i_poly[1] ^ i_poly[8];
			shift_poly_18_10[6] = i_poly[8] ^ i_poly[2] ^ i_poly[9];
			shift_poly_18_10[7] = i_poly[9] ^ i_poly[3];
			shift_poly_18_10[8] = i_poly[0] ^ i_poly[4];
			shift_poly_18_10[9] = i_poly[1] ^ i_poly[5];
			shift_poly_18_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_19_10;
		input [9:0] i_poly;
		begin
			shift_poly_19_10[0] = i_poly[1] ^ i_poly[5];
			shift_poly_19_10[1] = i_poly[2] ^ i_poly[6];
			shift_poly_19_10[2] = i_poly[3] ^ i_poly[0] ^ i_poly[7];
			shift_poly_19_10[3] = i_poly[4] ^ i_poly[5] ^ i_poly[8];
			shift_poly_19_10[4] = i_poly[5] ^ i_poly[9] ^ i_poly[6];
			shift_poly_19_10[5] = i_poly[6] ^ i_poly[0] ^ i_poly[7];
			shift_poly_19_10[6] = i_poly[7] ^ i_poly[1] ^ i_poly[8];
			shift_poly_19_10[7] = i_poly[8] ^ i_poly[2] ^ i_poly[9];
			shift_poly_19_10[8] = i_poly[9] ^ i_poly[3];
			shift_poly_19_10[9] = i_poly[0] ^ i_poly[4];
			shift_poly_19_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_20_10;
		input [9:0] i_poly;
		begin
			shift_poly_20_10[0] = i_poly[0] ^ i_poly[4];
			shift_poly_20_10[1] = i_poly[1] ^ i_poly[5];
			shift_poly_20_10[2] = i_poly[2] ^ i_poly[6];
			shift_poly_20_10[3] = i_poly[3] ^ i_poly[4] ^ i_poly[7];
			shift_poly_20_10[4] = i_poly[4] ^ i_poly[5] ^ i_poly[8];
			shift_poly_20_10[5] = i_poly[5] ^ i_poly[9] ^ i_poly[6];
			shift_poly_20_10[6] = i_poly[6] ^ i_poly[0] ^ i_poly[7];
			shift_poly_20_10[7] = i_poly[7] ^ i_poly[1] ^ i_poly[8];
			shift_poly_20_10[8] = i_poly[8] ^ i_poly[2] ^ i_poly[9];
			shift_poly_20_10[9] = i_poly[9] ^ i_poly[3];
			shift_poly_20_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_21_10;
		input [9:0] i_poly;
		begin
			shift_poly_21_10[0] = i_poly[9] ^ i_poly[3];
			shift_poly_21_10[1] = i_poly[0] ^ i_poly[4];
			shift_poly_21_10[2] = i_poly[1] ^ i_poly[5];
			shift_poly_21_10[3] = i_poly[2] ^ i_poly[6] ^ i_poly[9] ^ i_poly[3];
			shift_poly_21_10[4] = i_poly[3] ^ i_poly[4] ^ i_poly[7];
			shift_poly_21_10[5] = i_poly[4] ^ i_poly[5] ^ i_poly[8];
			shift_poly_21_10[6] = i_poly[5] ^ i_poly[9] ^ i_poly[6];
			shift_poly_21_10[7] = i_poly[6] ^ i_poly[0] ^ i_poly[7];
			shift_poly_21_10[8] = i_poly[7] ^ i_poly[1] ^ i_poly[8];
			shift_poly_21_10[9] = i_poly[8] ^ i_poly[2] ^ i_poly[9];
			shift_poly_21_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_22_10;
		input [9:0] i_poly;
		begin
			shift_poly_22_10[0] = i_poly[8] ^ i_poly[2] ^ i_poly[9];
			shift_poly_22_10[1] = i_poly[9] ^ i_poly[3];
			shift_poly_22_10[2] = i_poly[0] ^ i_poly[4];
			shift_poly_22_10[3] = i_poly[1] ^ i_poly[5] ^ i_poly[8] ^ i_poly[2] ^ i_poly[9];
			shift_poly_22_10[4] = i_poly[2] ^ i_poly[6] ^ i_poly[9] ^ i_poly[3];
			shift_poly_22_10[5] = i_poly[3] ^ i_poly[4] ^ i_poly[7];
			shift_poly_22_10[6] = i_poly[4] ^ i_poly[5] ^ i_poly[8];
			shift_poly_22_10[7] = i_poly[5] ^ i_poly[9] ^ i_poly[6];
			shift_poly_22_10[8] = i_poly[6] ^ i_poly[0] ^ i_poly[7];
			shift_poly_22_10[9] = i_poly[7] ^ i_poly[1] ^ i_poly[8];
			shift_poly_22_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_23_10;
		input [9:0] i_poly;
		begin
			shift_poly_23_10[0] = i_poly[7] ^ i_poly[1] ^ i_poly[8];
			shift_poly_23_10[1] = i_poly[8] ^ i_poly[2] ^ i_poly[9];
			shift_poly_23_10[2] = i_poly[9] ^ i_poly[3];
			shift_poly_23_10[3] = i_poly[0] ^ i_poly[4] ^ i_poly[7] ^ i_poly[1] ^ i_poly[8];
			shift_poly_23_10[4] = i_poly[1] ^ i_poly[5] ^ i_poly[8] ^ i_poly[2] ^ i_poly[9];
			shift_poly_23_10[5] = i_poly[2] ^ i_poly[6] ^ i_poly[9] ^ i_poly[3];
			shift_poly_23_10[6] = i_poly[3] ^ i_poly[4] ^ i_poly[7];
			shift_poly_23_10[7] = i_poly[4] ^ i_poly[5] ^ i_poly[8];
			shift_poly_23_10[8] = i_poly[5] ^ i_poly[9] ^ i_poly[6];
			shift_poly_23_10[9] = i_poly[6] ^ i_poly[0] ^ i_poly[7];
			shift_poly_23_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_24_10;
		input [9:0] i_poly;
		begin
			shift_poly_24_10[0] = i_poly[6] ^ i_poly[0] ^ i_poly[7];
			shift_poly_24_10[1] = i_poly[7] ^ i_poly[1] ^ i_poly[8];
			shift_poly_24_10[2] = i_poly[8] ^ i_poly[2] ^ i_poly[9];
			shift_poly_24_10[3] = i_poly[9] ^ i_poly[3] ^ i_poly[6] ^ i_poly[0] ^ i_poly[7];
			shift_poly_24_10[4] = i_poly[0] ^ i_poly[4] ^ i_poly[7] ^ i_poly[1] ^ i_poly[8];
			shift_poly_24_10[5] = i_poly[1] ^ i_poly[5] ^ i_poly[8] ^ i_poly[2] ^ i_poly[9];
			shift_poly_24_10[6] = i_poly[2] ^ i_poly[6] ^ i_poly[9] ^ i_poly[3];
			shift_poly_24_10[7] = i_poly[3] ^ i_poly[4] ^ i_poly[7];
			shift_poly_24_10[8] = i_poly[4] ^ i_poly[5] ^ i_poly[8];
			shift_poly_24_10[9] = i_poly[5] ^ i_poly[9] ^ i_poly[6];
			shift_poly_24_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_25_10;
		input [9:0] i_poly;
		begin
			shift_poly_25_10[0] = i_poly[5] ^ i_poly[9] ^ i_poly[6];
			shift_poly_25_10[1] = i_poly[6] ^ i_poly[0] ^ i_poly[7];
			shift_poly_25_10[2] = i_poly[7] ^ i_poly[1] ^ i_poly[8];
			shift_poly_25_10[3] = i_poly[8] ^ i_poly[2] ^ i_poly[5] ^ i_poly[6];
			shift_poly_25_10[4] = i_poly[9] ^ i_poly[3] ^ i_poly[6] ^ i_poly[0] ^ i_poly[7];
			shift_poly_25_10[5] = i_poly[0] ^ i_poly[4] ^ i_poly[7] ^ i_poly[1] ^ i_poly[8];
			shift_poly_25_10[6] = i_poly[1] ^ i_poly[5] ^ i_poly[8] ^ i_poly[2] ^ i_poly[9];
			shift_poly_25_10[7] = i_poly[2] ^ i_poly[6] ^ i_poly[9] ^ i_poly[3];
			shift_poly_25_10[8] = i_poly[3] ^ i_poly[4] ^ i_poly[7];
			shift_poly_25_10[9] = i_poly[4] ^ i_poly[5] ^ i_poly[8];
			shift_poly_25_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_26_10;
		input [9:0] i_poly;
		begin
			shift_poly_26_10[0] = i_poly[4] ^ i_poly[5] ^ i_poly[8];
			shift_poly_26_10[1] = i_poly[5] ^ i_poly[9] ^ i_poly[6];
			shift_poly_26_10[2] = i_poly[6] ^ i_poly[0] ^ i_poly[7];
			shift_poly_26_10[3] = i_poly[7] ^ i_poly[1] ^ i_poly[4] ^ i_poly[5];
			shift_poly_26_10[4] = i_poly[8] ^ i_poly[2] ^ i_poly[5] ^ i_poly[6];
			shift_poly_26_10[5] = i_poly[9] ^ i_poly[3] ^ i_poly[6] ^ i_poly[0] ^ i_poly[7];
			shift_poly_26_10[6] = i_poly[0] ^ i_poly[4] ^ i_poly[7] ^ i_poly[1] ^ i_poly[8];
			shift_poly_26_10[7] = i_poly[1] ^ i_poly[5] ^ i_poly[8] ^ i_poly[2] ^ i_poly[9];
			shift_poly_26_10[8] = i_poly[2] ^ i_poly[6] ^ i_poly[9] ^ i_poly[3];
			shift_poly_26_10[9] = i_poly[3] ^ i_poly[4] ^ i_poly[7];
			shift_poly_26_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_27_10;
		input [9:0] i_poly;
		begin
			shift_poly_27_10[0] = i_poly[3] ^ i_poly[4] ^ i_poly[7];
			shift_poly_27_10[1] = i_poly[4] ^ i_poly[5] ^ i_poly[8];
			shift_poly_27_10[2] = i_poly[5] ^ i_poly[9] ^ i_poly[6];
			shift_poly_27_10[3] = i_poly[6] ^ i_poly[0] ^ i_poly[3] ^ i_poly[4];
			shift_poly_27_10[4] = i_poly[7] ^ i_poly[1] ^ i_poly[4] ^ i_poly[5];
			shift_poly_27_10[5] = i_poly[8] ^ i_poly[2] ^ i_poly[5] ^ i_poly[6];
			shift_poly_27_10[6] = i_poly[9] ^ i_poly[3] ^ i_poly[6] ^ i_poly[0] ^ i_poly[7];
			shift_poly_27_10[7] = i_poly[0] ^ i_poly[4] ^ i_poly[7] ^ i_poly[1] ^ i_poly[8];
			shift_poly_27_10[8] = i_poly[1] ^ i_poly[5] ^ i_poly[8] ^ i_poly[2] ^ i_poly[9];
			shift_poly_27_10[9] = i_poly[2] ^ i_poly[6] ^ i_poly[9] ^ i_poly[3];
			shift_poly_27_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_28_10;
		input [9:0] i_poly;
		begin
			shift_poly_28_10[0] = i_poly[2] ^ i_poly[6] ^ i_poly[9] ^ i_poly[3];
			shift_poly_28_10[1] = i_poly[3] ^ i_poly[4] ^ i_poly[7];
			shift_poly_28_10[2] = i_poly[4] ^ i_poly[5] ^ i_poly[8];
			shift_poly_28_10[3] = i_poly[5] ^ i_poly[2] ^ i_poly[3];
			shift_poly_28_10[4] = i_poly[6] ^ i_poly[0] ^ i_poly[3] ^ i_poly[4];
			shift_poly_28_10[5] = i_poly[7] ^ i_poly[1] ^ i_poly[4] ^ i_poly[5];
			shift_poly_28_10[6] = i_poly[8] ^ i_poly[2] ^ i_poly[5] ^ i_poly[6];
			shift_poly_28_10[7] = i_poly[9] ^ i_poly[3] ^ i_poly[6] ^ i_poly[0] ^ i_poly[7];
			shift_poly_28_10[8] = i_poly[0] ^ i_poly[4] ^ i_poly[7] ^ i_poly[1] ^ i_poly[8];
			shift_poly_28_10[9] = i_poly[1] ^ i_poly[5] ^ i_poly[8] ^ i_poly[2] ^ i_poly[9];
			shift_poly_28_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_29_10;
		input [9:0] i_poly;
		begin
			shift_poly_29_10[0] = i_poly[1] ^ i_poly[5] ^ i_poly[8] ^ i_poly[2] ^ i_poly[9];
			shift_poly_29_10[1] = i_poly[2] ^ i_poly[6] ^ i_poly[9] ^ i_poly[3];
			shift_poly_29_10[2] = i_poly[3] ^ i_poly[4] ^ i_poly[7];
			shift_poly_29_10[3] = i_poly[4] ^ i_poly[1] ^ i_poly[2] ^ i_poly[9];
			shift_poly_29_10[4] = i_poly[5] ^ i_poly[2] ^ i_poly[3];
			shift_poly_29_10[5] = i_poly[6] ^ i_poly[0] ^ i_poly[3] ^ i_poly[4];
			shift_poly_29_10[6] = i_poly[7] ^ i_poly[1] ^ i_poly[4] ^ i_poly[5];
			shift_poly_29_10[7] = i_poly[8] ^ i_poly[2] ^ i_poly[5] ^ i_poly[6];
			shift_poly_29_10[8] = i_poly[9] ^ i_poly[3] ^ i_poly[6] ^ i_poly[0] ^ i_poly[7];
			shift_poly_29_10[9] = i_poly[0] ^ i_poly[4] ^ i_poly[7] ^ i_poly[1] ^ i_poly[8];
			shift_poly_29_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_30_10;
		input [9:0] i_poly;
		begin
			shift_poly_30_10[0] = i_poly[0] ^ i_poly[4] ^ i_poly[7] ^ i_poly[1] ^ i_poly[8];
			shift_poly_30_10[1] = i_poly[1] ^ i_poly[5] ^ i_poly[8] ^ i_poly[2] ^ i_poly[9];
			shift_poly_30_10[2] = i_poly[2] ^ i_poly[6] ^ i_poly[9] ^ i_poly[3];
			shift_poly_30_10[3] = i_poly[3] ^ i_poly[0] ^ i_poly[1] ^ i_poly[8];
			shift_poly_30_10[4] = i_poly[4] ^ i_poly[1] ^ i_poly[2] ^ i_poly[9];
			shift_poly_30_10[5] = i_poly[5] ^ i_poly[2] ^ i_poly[3];
			shift_poly_30_10[6] = i_poly[6] ^ i_poly[0] ^ i_poly[3] ^ i_poly[4];
			shift_poly_30_10[7] = i_poly[7] ^ i_poly[1] ^ i_poly[4] ^ i_poly[5];
			shift_poly_30_10[8] = i_poly[8] ^ i_poly[2] ^ i_poly[5] ^ i_poly[6];
			shift_poly_30_10[9] = i_poly[9] ^ i_poly[3] ^ i_poly[6] ^ i_poly[0] ^ i_poly[7];
			shift_poly_30_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_31_10;
		input [9:0] i_poly;
		begin
			shift_poly_31_10[0] = i_poly[9] ^ i_poly[3] ^ i_poly[6] ^ i_poly[0] ^ i_poly[7];
			shift_poly_31_10[1] = i_poly[0] ^ i_poly[4] ^ i_poly[7] ^ i_poly[1] ^ i_poly[8];
			shift_poly_31_10[2] = i_poly[1] ^ i_poly[5] ^ i_poly[8] ^ i_poly[2] ^ i_poly[9];
			shift_poly_31_10[3] = i_poly[2] ^ i_poly[0] ^ i_poly[7];
			shift_poly_31_10[4] = i_poly[3] ^ i_poly[0] ^ i_poly[1] ^ i_poly[8];
			shift_poly_31_10[5] = i_poly[4] ^ i_poly[1] ^ i_poly[2] ^ i_poly[9];
			shift_poly_31_10[6] = i_poly[5] ^ i_poly[2] ^ i_poly[3];
			shift_poly_31_10[7] = i_poly[6] ^ i_poly[0] ^ i_poly[3] ^ i_poly[4];
			shift_poly_31_10[8] = i_poly[7] ^ i_poly[1] ^ i_poly[4] ^ i_poly[5];
			shift_poly_31_10[9] = i_poly[8] ^ i_poly[2] ^ i_poly[5] ^ i_poly[6];
			shift_poly_31_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_32_10;
		input [9:0] i_poly;
		begin
			shift_poly_32_10[0] = i_poly[8] ^ i_poly[2] ^ i_poly[5] ^ i_poly[6];
			shift_poly_32_10[1] = i_poly[9] ^ i_poly[3] ^ i_poly[6] ^ i_poly[0] ^ i_poly[7];
			shift_poly_32_10[2] = i_poly[0] ^ i_poly[4] ^ i_poly[7] ^ i_poly[1] ^ i_poly[8];
			shift_poly_32_10[3] = i_poly[1] ^ i_poly[9] ^ i_poly[6];
			shift_poly_32_10[4] = i_poly[2] ^ i_poly[0] ^ i_poly[7];
			shift_poly_32_10[5] = i_poly[3] ^ i_poly[0] ^ i_poly[1] ^ i_poly[8];
			shift_poly_32_10[6] = i_poly[4] ^ i_poly[1] ^ i_poly[2] ^ i_poly[9];
			shift_poly_32_10[7] = i_poly[5] ^ i_poly[2] ^ i_poly[3];
			shift_poly_32_10[8] = i_poly[6] ^ i_poly[0] ^ i_poly[3] ^ i_poly[4];
			shift_poly_32_10[9] = i_poly[7] ^ i_poly[1] ^ i_poly[4] ^ i_poly[5];
			shift_poly_32_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_33_10;
		input [9:0] i_poly;
		begin
			shift_poly_33_10[0] = i_poly[7] ^ i_poly[1] ^ i_poly[4] ^ i_poly[5];
			shift_poly_33_10[1] = i_poly[8] ^ i_poly[2] ^ i_poly[5] ^ i_poly[6];
			shift_poly_33_10[2] = i_poly[9] ^ i_poly[3] ^ i_poly[6] ^ i_poly[0] ^ i_poly[7];
			shift_poly_33_10[3] = i_poly[0] ^ i_poly[8] ^ i_poly[5];
			shift_poly_33_10[4] = i_poly[1] ^ i_poly[9] ^ i_poly[6];
			shift_poly_33_10[5] = i_poly[2] ^ i_poly[0] ^ i_poly[7];
			shift_poly_33_10[6] = i_poly[3] ^ i_poly[0] ^ i_poly[1] ^ i_poly[8];
			shift_poly_33_10[7] = i_poly[4] ^ i_poly[1] ^ i_poly[2] ^ i_poly[9];
			shift_poly_33_10[8] = i_poly[5] ^ i_poly[2] ^ i_poly[3];
			shift_poly_33_10[9] = i_poly[6] ^ i_poly[0] ^ i_poly[3] ^ i_poly[4];
			shift_poly_33_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_34_10;
		input [9:0] i_poly;
		begin
			shift_poly_34_10[0] = i_poly[6] ^ i_poly[0] ^ i_poly[3] ^ i_poly[4];
			shift_poly_34_10[1] = i_poly[7] ^ i_poly[1] ^ i_poly[4] ^ i_poly[5];
			shift_poly_34_10[2] = i_poly[8] ^ i_poly[2] ^ i_poly[5] ^ i_poly[6];
			shift_poly_34_10[3] = i_poly[9] ^ i_poly[7] ^ i_poly[4];
			shift_poly_34_10[4] = i_poly[0] ^ i_poly[8] ^ i_poly[5];
			shift_poly_34_10[5] = i_poly[1] ^ i_poly[9] ^ i_poly[6];
			shift_poly_34_10[6] = i_poly[2] ^ i_poly[0] ^ i_poly[7];
			shift_poly_34_10[7] = i_poly[3] ^ i_poly[0] ^ i_poly[1] ^ i_poly[8];
			shift_poly_34_10[8] = i_poly[4] ^ i_poly[1] ^ i_poly[2] ^ i_poly[9];
			shift_poly_34_10[9] = i_poly[5] ^ i_poly[2] ^ i_poly[3];
			shift_poly_34_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_35_10;
		input [9:0] i_poly;
		begin
			shift_poly_35_10[0] = i_poly[5] ^ i_poly[2] ^ i_poly[3];
			shift_poly_35_10[1] = i_poly[6] ^ i_poly[0] ^ i_poly[3] ^ i_poly[4];
			shift_poly_35_10[2] = i_poly[7] ^ i_poly[1] ^ i_poly[4] ^ i_poly[5];
			shift_poly_35_10[3] = i_poly[8] ^ i_poly[6] ^ i_poly[3];
			shift_poly_35_10[4] = i_poly[9] ^ i_poly[7] ^ i_poly[4];
			shift_poly_35_10[5] = i_poly[0] ^ i_poly[8] ^ i_poly[5];
			shift_poly_35_10[6] = i_poly[1] ^ i_poly[9] ^ i_poly[6];
			shift_poly_35_10[7] = i_poly[2] ^ i_poly[0] ^ i_poly[7];
			shift_poly_35_10[8] = i_poly[3] ^ i_poly[0] ^ i_poly[1] ^ i_poly[8];
			shift_poly_35_10[9] = i_poly[4] ^ i_poly[1] ^ i_poly[2] ^ i_poly[9];
			shift_poly_35_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_36_10;
		input [9:0] i_poly;
		begin
			shift_poly_36_10[0] = i_poly[4] ^ i_poly[1] ^ i_poly[2] ^ i_poly[9];
			shift_poly_36_10[1] = i_poly[5] ^ i_poly[2] ^ i_poly[3];
			shift_poly_36_10[2] = i_poly[6] ^ i_poly[0] ^ i_poly[3] ^ i_poly[4];
			shift_poly_36_10[3] = i_poly[7] ^ i_poly[5] ^ i_poly[9] ^ i_poly[2];
			shift_poly_36_10[4] = i_poly[8] ^ i_poly[6] ^ i_poly[3];
			shift_poly_36_10[5] = i_poly[9] ^ i_poly[7] ^ i_poly[4];
			shift_poly_36_10[6] = i_poly[0] ^ i_poly[8] ^ i_poly[5];
			shift_poly_36_10[7] = i_poly[1] ^ i_poly[9] ^ i_poly[6];
			shift_poly_36_10[8] = i_poly[2] ^ i_poly[0] ^ i_poly[7];
			shift_poly_36_10[9] = i_poly[3] ^ i_poly[0] ^ i_poly[1] ^ i_poly[8];
			shift_poly_36_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_37_10;
		input [9:0] i_poly;
		begin
			shift_poly_37_10[0] = i_poly[3] ^ i_poly[0] ^ i_poly[1] ^ i_poly[8];
			shift_poly_37_10[1] = i_poly[4] ^ i_poly[1] ^ i_poly[2] ^ i_poly[9];
			shift_poly_37_10[2] = i_poly[5] ^ i_poly[2] ^ i_poly[3];
			shift_poly_37_10[3] = i_poly[6] ^ i_poly[4] ^ i_poly[1] ^ i_poly[8];
			shift_poly_37_10[4] = i_poly[7] ^ i_poly[5] ^ i_poly[2] ^ i_poly[9];
			shift_poly_37_10[5] = i_poly[8] ^ i_poly[6] ^ i_poly[3];
			shift_poly_37_10[6] = i_poly[9] ^ i_poly[7] ^ i_poly[4];
			shift_poly_37_10[7] = i_poly[0] ^ i_poly[8] ^ i_poly[5];
			shift_poly_37_10[8] = i_poly[1] ^ i_poly[9] ^ i_poly[6];
			shift_poly_37_10[9] = i_poly[2] ^ i_poly[0] ^ i_poly[7];
			shift_poly_37_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_38_10;
		input [9:0] i_poly;
		begin
			shift_poly_38_10[0] = i_poly[2] ^ i_poly[0] ^ i_poly[7];
			shift_poly_38_10[1] = i_poly[3] ^ i_poly[0] ^ i_poly[1] ^ i_poly[8];
			shift_poly_38_10[2] = i_poly[4] ^ i_poly[1] ^ i_poly[2] ^ i_poly[9];
			shift_poly_38_10[3] = i_poly[5] ^ i_poly[3] ^ i_poly[0] ^ i_poly[7];
			shift_poly_38_10[4] = i_poly[6] ^ i_poly[4] ^ i_poly[1] ^ i_poly[8];
			shift_poly_38_10[5] = i_poly[7] ^ i_poly[5] ^ i_poly[2] ^ i_poly[9];
			shift_poly_38_10[6] = i_poly[8] ^ i_poly[6] ^ i_poly[3];
			shift_poly_38_10[7] = i_poly[9] ^ i_poly[7] ^ i_poly[4];
			shift_poly_38_10[8] = i_poly[0] ^ i_poly[8] ^ i_poly[5];
			shift_poly_38_10[9] = i_poly[1] ^ i_poly[9] ^ i_poly[6];
			shift_poly_38_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_39_10;
		input [9:0] i_poly;
		begin
			shift_poly_39_10[0] = i_poly[1] ^ i_poly[9] ^ i_poly[6];
			shift_poly_39_10[1] = i_poly[2] ^ i_poly[0] ^ i_poly[7];
			shift_poly_39_10[2] = i_poly[3] ^ i_poly[0] ^ i_poly[1] ^ i_poly[8];
			shift_poly_39_10[3] = i_poly[4] ^ i_poly[2] ^ i_poly[6];
			shift_poly_39_10[4] = i_poly[5] ^ i_poly[3] ^ i_poly[0] ^ i_poly[7];
			shift_poly_39_10[5] = i_poly[6] ^ i_poly[4] ^ i_poly[1] ^ i_poly[8];
			shift_poly_39_10[6] = i_poly[7] ^ i_poly[5] ^ i_poly[2] ^ i_poly[9];
			shift_poly_39_10[7] = i_poly[8] ^ i_poly[6] ^ i_poly[3];
			shift_poly_39_10[8] = i_poly[9] ^ i_poly[7] ^ i_poly[4];
			shift_poly_39_10[9] = i_poly[0] ^ i_poly[8] ^ i_poly[5];
			shift_poly_39_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_40_10;
		input [9:0] i_poly;
		begin
			shift_poly_40_10[0] = i_poly[0] ^ i_poly[8] ^ i_poly[5];
			shift_poly_40_10[1] = i_poly[1] ^ i_poly[9] ^ i_poly[6];
			shift_poly_40_10[2] = i_poly[2] ^ i_poly[0] ^ i_poly[7];
			shift_poly_40_10[3] = i_poly[3] ^ i_poly[1] ^ i_poly[5];
			shift_poly_40_10[4] = i_poly[4] ^ i_poly[2] ^ i_poly[6];
			shift_poly_40_10[5] = i_poly[5] ^ i_poly[3] ^ i_poly[0] ^ i_poly[7];
			shift_poly_40_10[6] = i_poly[6] ^ i_poly[4] ^ i_poly[1] ^ i_poly[8];
			shift_poly_40_10[7] = i_poly[7] ^ i_poly[5] ^ i_poly[2] ^ i_poly[9];
			shift_poly_40_10[8] = i_poly[8] ^ i_poly[6] ^ i_poly[3];
			shift_poly_40_10[9] = i_poly[9] ^ i_poly[7] ^ i_poly[4];
			shift_poly_40_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_41_10;
		input [9:0] i_poly;
		begin
			shift_poly_41_10[0] = i_poly[9] ^ i_poly[7] ^ i_poly[4];
			shift_poly_41_10[1] = i_poly[0] ^ i_poly[8] ^ i_poly[5];
			shift_poly_41_10[2] = i_poly[1] ^ i_poly[9] ^ i_poly[6];
			shift_poly_41_10[3] = i_poly[2] ^ i_poly[0] ^ i_poly[4] ^ i_poly[9];
			shift_poly_41_10[4] = i_poly[3] ^ i_poly[1] ^ i_poly[5];
			shift_poly_41_10[5] = i_poly[4] ^ i_poly[2] ^ i_poly[6];
			shift_poly_41_10[6] = i_poly[5] ^ i_poly[3] ^ i_poly[0] ^ i_poly[7];
			shift_poly_41_10[7] = i_poly[6] ^ i_poly[4] ^ i_poly[1] ^ i_poly[8];
			shift_poly_41_10[8] = i_poly[7] ^ i_poly[5] ^ i_poly[2] ^ i_poly[9];
			shift_poly_41_10[9] = i_poly[8] ^ i_poly[6] ^ i_poly[3];
			shift_poly_41_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_42_10;
		input [9:0] i_poly;
		begin
			shift_poly_42_10[0] = i_poly[8] ^ i_poly[6] ^ i_poly[3];
			shift_poly_42_10[1] = i_poly[9] ^ i_poly[7] ^ i_poly[4];
			shift_poly_42_10[2] = i_poly[0] ^ i_poly[8] ^ i_poly[5];
			shift_poly_42_10[3] = i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[8];
			shift_poly_42_10[4] = i_poly[2] ^ i_poly[0] ^ i_poly[4] ^ i_poly[9];
			shift_poly_42_10[5] = i_poly[3] ^ i_poly[1] ^ i_poly[5];
			shift_poly_42_10[6] = i_poly[4] ^ i_poly[2] ^ i_poly[6];
			shift_poly_42_10[7] = i_poly[5] ^ i_poly[3] ^ i_poly[0] ^ i_poly[7];
			shift_poly_42_10[8] = i_poly[6] ^ i_poly[4] ^ i_poly[1] ^ i_poly[8];
			shift_poly_42_10[9] = i_poly[7] ^ i_poly[5] ^ i_poly[2] ^ i_poly[9];
			shift_poly_42_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_43_10;
		input [9:0] i_poly;
		begin
			shift_poly_43_10[0] = i_poly[7] ^ i_poly[5] ^ i_poly[2] ^ i_poly[9];
			shift_poly_43_10[1] = i_poly[8] ^ i_poly[6] ^ i_poly[3];
			shift_poly_43_10[2] = i_poly[9] ^ i_poly[7] ^ i_poly[4];
			shift_poly_43_10[3] = i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_43_10[4] = i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[8];
			shift_poly_43_10[5] = i_poly[2] ^ i_poly[0] ^ i_poly[4] ^ i_poly[9];
			shift_poly_43_10[6] = i_poly[3] ^ i_poly[1] ^ i_poly[5];
			shift_poly_43_10[7] = i_poly[4] ^ i_poly[2] ^ i_poly[6];
			shift_poly_43_10[8] = i_poly[5] ^ i_poly[3] ^ i_poly[0] ^ i_poly[7];
			shift_poly_43_10[9] = i_poly[6] ^ i_poly[4] ^ i_poly[1] ^ i_poly[8];
			shift_poly_43_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_44_10;
		input [9:0] i_poly;
		begin
			shift_poly_44_10[0] = i_poly[6] ^ i_poly[4] ^ i_poly[1] ^ i_poly[8];
			shift_poly_44_10[1] = i_poly[7] ^ i_poly[5] ^ i_poly[2] ^ i_poly[9];
			shift_poly_44_10[2] = i_poly[8] ^ i_poly[6] ^ i_poly[3];
			shift_poly_44_10[3] = i_poly[9] ^ i_poly[7] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_44_10[4] = i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_44_10[5] = i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[8];
			shift_poly_44_10[6] = i_poly[2] ^ i_poly[0] ^ i_poly[4] ^ i_poly[9];
			shift_poly_44_10[7] = i_poly[3] ^ i_poly[1] ^ i_poly[5];
			shift_poly_44_10[8] = i_poly[4] ^ i_poly[2] ^ i_poly[6];
			shift_poly_44_10[9] = i_poly[5] ^ i_poly[3] ^ i_poly[0] ^ i_poly[7];
			shift_poly_44_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_45_10;
		input [9:0] i_poly;
		begin
			shift_poly_45_10[0] = i_poly[5] ^ i_poly[3] ^ i_poly[0] ^ i_poly[7];
			shift_poly_45_10[1] = i_poly[6] ^ i_poly[4] ^ i_poly[1] ^ i_poly[8];
			shift_poly_45_10[2] = i_poly[7] ^ i_poly[5] ^ i_poly[2] ^ i_poly[9];
			shift_poly_45_10[3] = i_poly[8] ^ i_poly[6] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_45_10[4] = i_poly[9] ^ i_poly[7] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_45_10[5] = i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_45_10[6] = i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[8];
			shift_poly_45_10[7] = i_poly[2] ^ i_poly[0] ^ i_poly[4] ^ i_poly[9];
			shift_poly_45_10[8] = i_poly[3] ^ i_poly[1] ^ i_poly[5];
			shift_poly_45_10[9] = i_poly[4] ^ i_poly[2] ^ i_poly[6];
			shift_poly_45_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_46_10;
		input [9:0] i_poly;
		begin
			shift_poly_46_10[0] = i_poly[4] ^ i_poly[2] ^ i_poly[6];
			shift_poly_46_10[1] = i_poly[5] ^ i_poly[3] ^ i_poly[0] ^ i_poly[7];
			shift_poly_46_10[2] = i_poly[6] ^ i_poly[4] ^ i_poly[1] ^ i_poly[8];
			shift_poly_46_10[3] = i_poly[7] ^ i_poly[5] ^ i_poly[9] ^ i_poly[4] ^ i_poly[6];
			shift_poly_46_10[4] = i_poly[8] ^ i_poly[6] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_46_10[5] = i_poly[9] ^ i_poly[7] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_46_10[6] = i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_46_10[7] = i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[8];
			shift_poly_46_10[8] = i_poly[2] ^ i_poly[0] ^ i_poly[4] ^ i_poly[9];
			shift_poly_46_10[9] = i_poly[3] ^ i_poly[1] ^ i_poly[5];
			shift_poly_46_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_47_10;
		input [9:0] i_poly;
		begin
			shift_poly_47_10[0] = i_poly[3] ^ i_poly[1] ^ i_poly[5];
			shift_poly_47_10[1] = i_poly[4] ^ i_poly[2] ^ i_poly[6];
			shift_poly_47_10[2] = i_poly[5] ^ i_poly[3] ^ i_poly[0] ^ i_poly[7];
			shift_poly_47_10[3] = i_poly[6] ^ i_poly[4] ^ i_poly[8] ^ i_poly[3] ^ i_poly[5];
			shift_poly_47_10[4] = i_poly[7] ^ i_poly[5] ^ i_poly[9] ^ i_poly[4] ^ i_poly[6];
			shift_poly_47_10[5] = i_poly[8] ^ i_poly[6] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_47_10[6] = i_poly[9] ^ i_poly[7] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_47_10[7] = i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_47_10[8] = i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[8];
			shift_poly_47_10[9] = i_poly[2] ^ i_poly[0] ^ i_poly[4] ^ i_poly[9];
			shift_poly_47_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_48_10;
		input [9:0] i_poly;
		begin
			shift_poly_48_10[0] = i_poly[2] ^ i_poly[0] ^ i_poly[4] ^ i_poly[9];
			shift_poly_48_10[1] = i_poly[3] ^ i_poly[1] ^ i_poly[5];
			shift_poly_48_10[2] = i_poly[4] ^ i_poly[2] ^ i_poly[6];
			shift_poly_48_10[3] = i_poly[5] ^ i_poly[3] ^ i_poly[7] ^ i_poly[2] ^ i_poly[4] ^ i_poly[9];
			shift_poly_48_10[4] = i_poly[6] ^ i_poly[4] ^ i_poly[8] ^ i_poly[3] ^ i_poly[5];
			shift_poly_48_10[5] = i_poly[7] ^ i_poly[5] ^ i_poly[9] ^ i_poly[4] ^ i_poly[6];
			shift_poly_48_10[6] = i_poly[8] ^ i_poly[6] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_48_10[7] = i_poly[9] ^ i_poly[7] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_48_10[8] = i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_48_10[9] = i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[8];
			shift_poly_48_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_49_10;
		input [9:0] i_poly;
		begin
			shift_poly_49_10[0] = i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[8];
			shift_poly_49_10[1] = i_poly[2] ^ i_poly[0] ^ i_poly[4] ^ i_poly[9];
			shift_poly_49_10[2] = i_poly[3] ^ i_poly[1] ^ i_poly[5];
			shift_poly_49_10[3] = i_poly[4] ^ i_poly[2] ^ i_poly[6] ^ i_poly[1] ^ i_poly[3] ^ i_poly[8] ^ i_poly[9];
			shift_poly_49_10[4] = i_poly[5] ^ i_poly[3] ^ i_poly[7] ^ i_poly[2] ^ i_poly[4] ^ i_poly[9];
			shift_poly_49_10[5] = i_poly[6] ^ i_poly[4] ^ i_poly[8] ^ i_poly[3] ^ i_poly[5];
			shift_poly_49_10[6] = i_poly[7] ^ i_poly[5] ^ i_poly[9] ^ i_poly[4] ^ i_poly[6];
			shift_poly_49_10[7] = i_poly[8] ^ i_poly[6] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_49_10[8] = i_poly[9] ^ i_poly[7] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_49_10[9] = i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_49_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_50_10;
		input [9:0] i_poly;
		begin
			shift_poly_50_10[0] = i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_50_10[1] = i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[8];
			shift_poly_50_10[2] = i_poly[2] ^ i_poly[0] ^ i_poly[4] ^ i_poly[9];
			shift_poly_50_10[3] = i_poly[3] ^ i_poly[1] ^ i_poly[5] ^ i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_50_10[4] = i_poly[4] ^ i_poly[2] ^ i_poly[6] ^ i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[8];
			shift_poly_50_10[5] = i_poly[5] ^ i_poly[3] ^ i_poly[7] ^ i_poly[2] ^ i_poly[4] ^ i_poly[9];
			shift_poly_50_10[6] = i_poly[6] ^ i_poly[4] ^ i_poly[8] ^ i_poly[3] ^ i_poly[5];
			shift_poly_50_10[7] = i_poly[7] ^ i_poly[5] ^ i_poly[9] ^ i_poly[4] ^ i_poly[6];
			shift_poly_50_10[8] = i_poly[8] ^ i_poly[6] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_50_10[9] = i_poly[9] ^ i_poly[7] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_50_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_51_10;
		input [9:0] i_poly;
		begin
			shift_poly_51_10[0] = i_poly[9] ^ i_poly[7] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_51_10[1] = i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_51_10[2] = i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[8];
			shift_poly_51_10[3] = i_poly[2] ^ i_poly[0] ^ i_poly[4] ^ i_poly[7] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_51_10[4] = i_poly[3] ^ i_poly[1] ^ i_poly[5] ^ i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_51_10[5] = i_poly[4] ^ i_poly[2] ^ i_poly[6] ^ i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[8];
			shift_poly_51_10[6] = i_poly[5] ^ i_poly[3] ^ i_poly[7] ^ i_poly[2] ^ i_poly[4] ^ i_poly[9];
			shift_poly_51_10[7] = i_poly[6] ^ i_poly[4] ^ i_poly[8] ^ i_poly[3] ^ i_poly[5];
			shift_poly_51_10[8] = i_poly[7] ^ i_poly[5] ^ i_poly[9] ^ i_poly[4] ^ i_poly[6];
			shift_poly_51_10[9] = i_poly[8] ^ i_poly[6] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_51_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_52_10;
		input [9:0] i_poly;
		begin
			shift_poly_52_10[0] = i_poly[8] ^ i_poly[6] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_52_10[1] = i_poly[9] ^ i_poly[7] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_52_10[2] = i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_52_10[3] = i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[6] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_52_10[4] = i_poly[2] ^ i_poly[0] ^ i_poly[4] ^ i_poly[7] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_52_10[5] = i_poly[3] ^ i_poly[1] ^ i_poly[5] ^ i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_52_10[6] = i_poly[4] ^ i_poly[2] ^ i_poly[6] ^ i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[8];
			shift_poly_52_10[7] = i_poly[5] ^ i_poly[3] ^ i_poly[7] ^ i_poly[2] ^ i_poly[4] ^ i_poly[9];
			shift_poly_52_10[8] = i_poly[6] ^ i_poly[4] ^ i_poly[8] ^ i_poly[3] ^ i_poly[5];
			shift_poly_52_10[9] = i_poly[7] ^ i_poly[5] ^ i_poly[9] ^ i_poly[4] ^ i_poly[6];
			shift_poly_52_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_53_10;
		input [9:0] i_poly;
		begin
			shift_poly_53_10[0] = i_poly[7] ^ i_poly[5] ^ i_poly[9] ^ i_poly[4] ^ i_poly[6];
			shift_poly_53_10[1] = i_poly[8] ^ i_poly[6] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_53_10[2] = i_poly[9] ^ i_poly[7] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_53_10[3] = i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[5] ^ i_poly[4] ^ i_poly[6];
			shift_poly_53_10[4] = i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[6] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_53_10[5] = i_poly[2] ^ i_poly[0] ^ i_poly[4] ^ i_poly[7] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_53_10[6] = i_poly[3] ^ i_poly[1] ^ i_poly[5] ^ i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_53_10[7] = i_poly[4] ^ i_poly[2] ^ i_poly[6] ^ i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[8];
			shift_poly_53_10[8] = i_poly[5] ^ i_poly[3] ^ i_poly[7] ^ i_poly[2] ^ i_poly[4] ^ i_poly[9];
			shift_poly_53_10[9] = i_poly[6] ^ i_poly[4] ^ i_poly[8] ^ i_poly[3] ^ i_poly[5];
			shift_poly_53_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_54_10;
		input [9:0] i_poly;
		begin
			shift_poly_54_10[0] = i_poly[6] ^ i_poly[4] ^ i_poly[8] ^ i_poly[3] ^ i_poly[5];
			shift_poly_54_10[1] = i_poly[7] ^ i_poly[5] ^ i_poly[9] ^ i_poly[4] ^ i_poly[6];
			shift_poly_54_10[2] = i_poly[8] ^ i_poly[6] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_54_10[3] = i_poly[9] ^ i_poly[7] ^ i_poly[1] ^ i_poly[4] ^ i_poly[3] ^ i_poly[5];
			shift_poly_54_10[4] = i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[5] ^ i_poly[4] ^ i_poly[6];
			shift_poly_54_10[5] = i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[6] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_54_10[6] = i_poly[2] ^ i_poly[0] ^ i_poly[4] ^ i_poly[7] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_54_10[7] = i_poly[3] ^ i_poly[1] ^ i_poly[5] ^ i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_54_10[8] = i_poly[4] ^ i_poly[2] ^ i_poly[6] ^ i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[8];
			shift_poly_54_10[9] = i_poly[5] ^ i_poly[3] ^ i_poly[7] ^ i_poly[2] ^ i_poly[4] ^ i_poly[9];
			shift_poly_54_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_55_10;
		input [9:0] i_poly;
		begin
			shift_poly_55_10[0] = i_poly[5] ^ i_poly[3] ^ i_poly[7] ^ i_poly[2] ^ i_poly[4] ^ i_poly[9];
			shift_poly_55_10[1] = i_poly[6] ^ i_poly[4] ^ i_poly[8] ^ i_poly[3] ^ i_poly[5];
			shift_poly_55_10[2] = i_poly[7] ^ i_poly[5] ^ i_poly[9] ^ i_poly[4] ^ i_poly[6];
			shift_poly_55_10[3] = i_poly[8] ^ i_poly[6] ^ i_poly[0] ^ i_poly[3] ^ i_poly[2] ^ i_poly[4] ^ i_poly[9];
			shift_poly_55_10[4] = i_poly[9] ^ i_poly[7] ^ i_poly[1] ^ i_poly[4] ^ i_poly[3] ^ i_poly[5];
			shift_poly_55_10[5] = i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[5] ^ i_poly[4] ^ i_poly[6];
			shift_poly_55_10[6] = i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[6] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_55_10[7] = i_poly[2] ^ i_poly[0] ^ i_poly[4] ^ i_poly[7] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_55_10[8] = i_poly[3] ^ i_poly[1] ^ i_poly[5] ^ i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_55_10[9] = i_poly[4] ^ i_poly[2] ^ i_poly[6] ^ i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[8];
			shift_poly_55_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_56_10;
		input [9:0] i_poly;
		begin
			shift_poly_56_10[0] = i_poly[4] ^ i_poly[2] ^ i_poly[6] ^ i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[8];
			shift_poly_56_10[1] = i_poly[5] ^ i_poly[3] ^ i_poly[7] ^ i_poly[2] ^ i_poly[4] ^ i_poly[9];
			shift_poly_56_10[2] = i_poly[6] ^ i_poly[4] ^ i_poly[8] ^ i_poly[3] ^ i_poly[5];
			shift_poly_56_10[3] = i_poly[7] ^ i_poly[5] ^ i_poly[2] ^ i_poly[1] ^ i_poly[3] ^ i_poly[8];
			shift_poly_56_10[4] = i_poly[8] ^ i_poly[6] ^ i_poly[0] ^ i_poly[3] ^ i_poly[2] ^ i_poly[4] ^ i_poly[9];
			shift_poly_56_10[5] = i_poly[9] ^ i_poly[7] ^ i_poly[1] ^ i_poly[4] ^ i_poly[3] ^ i_poly[5];
			shift_poly_56_10[6] = i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[5] ^ i_poly[4] ^ i_poly[6];
			shift_poly_56_10[7] = i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[6] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_56_10[8] = i_poly[2] ^ i_poly[0] ^ i_poly[4] ^ i_poly[7] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_56_10[9] = i_poly[3] ^ i_poly[1] ^ i_poly[5] ^ i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_56_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_57_10;
		input [9:0] i_poly;
		begin
			shift_poly_57_10[0] = i_poly[3] ^ i_poly[1] ^ i_poly[5] ^ i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_57_10[1] = i_poly[4] ^ i_poly[2] ^ i_poly[6] ^ i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[8];
			shift_poly_57_10[2] = i_poly[5] ^ i_poly[3] ^ i_poly[7] ^ i_poly[2] ^ i_poly[4] ^ i_poly[9];
			shift_poly_57_10[3] = i_poly[6] ^ i_poly[4] ^ i_poly[1] ^ i_poly[0] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_57_10[4] = i_poly[7] ^ i_poly[5] ^ i_poly[2] ^ i_poly[1] ^ i_poly[3] ^ i_poly[8];
			shift_poly_57_10[5] = i_poly[8] ^ i_poly[6] ^ i_poly[0] ^ i_poly[3] ^ i_poly[2] ^ i_poly[4] ^ i_poly[9];
			shift_poly_57_10[6] = i_poly[9] ^ i_poly[7] ^ i_poly[1] ^ i_poly[4] ^ i_poly[3] ^ i_poly[5];
			shift_poly_57_10[7] = i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[5] ^ i_poly[4] ^ i_poly[6];
			shift_poly_57_10[8] = i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[6] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_57_10[9] = i_poly[2] ^ i_poly[0] ^ i_poly[4] ^ i_poly[7] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_57_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_58_10;
		input [9:0] i_poly;
		begin
			shift_poly_58_10[0] = i_poly[2] ^ i_poly[0] ^ i_poly[4] ^ i_poly[7] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_58_10[1] = i_poly[3] ^ i_poly[1] ^ i_poly[5] ^ i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_58_10[2] = i_poly[4] ^ i_poly[2] ^ i_poly[6] ^ i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[8];
			shift_poly_58_10[3] = i_poly[5] ^ i_poly[3] ^ i_poly[9] ^ i_poly[0] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_58_10[4] = i_poly[6] ^ i_poly[4] ^ i_poly[1] ^ i_poly[0] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_58_10[5] = i_poly[7] ^ i_poly[5] ^ i_poly[2] ^ i_poly[1] ^ i_poly[3] ^ i_poly[8];
			shift_poly_58_10[6] = i_poly[8] ^ i_poly[6] ^ i_poly[0] ^ i_poly[3] ^ i_poly[2] ^ i_poly[4] ^ i_poly[9];
			shift_poly_58_10[7] = i_poly[9] ^ i_poly[7] ^ i_poly[1] ^ i_poly[4] ^ i_poly[3] ^ i_poly[5];
			shift_poly_58_10[8] = i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[5] ^ i_poly[4] ^ i_poly[6];
			shift_poly_58_10[9] = i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[6] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_58_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_59_10;
		input [9:0] i_poly;
		begin
			shift_poly_59_10[0] = i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[6] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_59_10[1] = i_poly[2] ^ i_poly[0] ^ i_poly[4] ^ i_poly[7] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_59_10[2] = i_poly[3] ^ i_poly[1] ^ i_poly[5] ^ i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_59_10[3] = i_poly[4] ^ i_poly[2] ^ i_poly[8] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_59_10[4] = i_poly[5] ^ i_poly[3] ^ i_poly[9] ^ i_poly[0] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_59_10[5] = i_poly[6] ^ i_poly[4] ^ i_poly[1] ^ i_poly[0] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_59_10[6] = i_poly[7] ^ i_poly[5] ^ i_poly[2] ^ i_poly[1] ^ i_poly[3] ^ i_poly[8];
			shift_poly_59_10[7] = i_poly[8] ^ i_poly[6] ^ i_poly[0] ^ i_poly[3] ^ i_poly[2] ^ i_poly[4] ^ i_poly[9];
			shift_poly_59_10[8] = i_poly[9] ^ i_poly[7] ^ i_poly[1] ^ i_poly[4] ^ i_poly[3] ^ i_poly[5];
			shift_poly_59_10[9] = i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[5] ^ i_poly[4] ^ i_poly[6];
			shift_poly_59_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_60_10;
		input [9:0] i_poly;
		begin
			shift_poly_60_10[0] = i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[5] ^ i_poly[4] ^ i_poly[6];
			shift_poly_60_10[1] = i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[6] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_60_10[2] = i_poly[2] ^ i_poly[0] ^ i_poly[4] ^ i_poly[7] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_60_10[3] = i_poly[3] ^ i_poly[1] ^ i_poly[7] ^ i_poly[9] ^ i_poly[4] ^ i_poly[6];
			shift_poly_60_10[4] = i_poly[4] ^ i_poly[2] ^ i_poly[8] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_60_10[5] = i_poly[5] ^ i_poly[3] ^ i_poly[9] ^ i_poly[0] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_60_10[6] = i_poly[6] ^ i_poly[4] ^ i_poly[1] ^ i_poly[0] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_60_10[7] = i_poly[7] ^ i_poly[5] ^ i_poly[2] ^ i_poly[1] ^ i_poly[3] ^ i_poly[8];
			shift_poly_60_10[8] = i_poly[8] ^ i_poly[6] ^ i_poly[0] ^ i_poly[3] ^ i_poly[2] ^ i_poly[4] ^ i_poly[9];
			shift_poly_60_10[9] = i_poly[9] ^ i_poly[7] ^ i_poly[1] ^ i_poly[4] ^ i_poly[3] ^ i_poly[5];
			shift_poly_60_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_61_10;
		input [9:0] i_poly;
		begin
			shift_poly_61_10[0] = i_poly[9] ^ i_poly[7] ^ i_poly[1] ^ i_poly[4] ^ i_poly[3] ^ i_poly[5];
			shift_poly_61_10[1] = i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[5] ^ i_poly[4] ^ i_poly[6];
			shift_poly_61_10[2] = i_poly[1] ^ i_poly[9] ^ i_poly[3] ^ i_poly[6] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_61_10[3] = i_poly[2] ^ i_poly[0] ^ i_poly[6] ^ i_poly[8] ^ i_poly[9] ^ i_poly[3] ^ i_poly[5];
			shift_poly_61_10[4] = i_poly[3] ^ i_poly[1] ^ i_poly[7] ^ i_poly[9] ^ i_poly[4] ^ i_poly[6];
			shift_poly_61_10[5] = i_poly[4] ^ i_poly[2] ^ i_poly[8] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_61_10[6] = i_poly[5] ^ i_poly[3] ^ i_poly[9] ^ i_poly[0] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_61_10[7] = i_poly[6] ^ i_poly[4] ^ i_poly[1] ^ i_poly[0] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_61_10[8] = i_poly[7] ^ i_poly[5] ^ i_poly[2] ^ i_poly[1] ^ i_poly[3] ^ i_poly[8];
			shift_poly_61_10[9] = i_poly[8] ^ i_poly[6] ^ i_poly[0] ^ i_poly[3] ^ i_poly[2] ^ i_poly[4] ^ i_poly[9];
			shift_poly_61_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_62_10;
		input [9:0] i_poly;
		begin
			shift_poly_62_10[0] = i_poly[8] ^ i_poly[6] ^ i_poly[0] ^ i_poly[3] ^ i_poly[2] ^ i_poly[4] ^ i_poly[9];
			shift_poly_62_10[1] = i_poly[9] ^ i_poly[7] ^ i_poly[1] ^ i_poly[4] ^ i_poly[3] ^ i_poly[5];
			shift_poly_62_10[2] = i_poly[0] ^ i_poly[8] ^ i_poly[2] ^ i_poly[5] ^ i_poly[4] ^ i_poly[6];
			shift_poly_62_10[3] = i_poly[1] ^ i_poly[5] ^ i_poly[7] ^ i_poly[8] ^ i_poly[2] ^ i_poly[4];
			shift_poly_62_10[4] = i_poly[2] ^ i_poly[0] ^ i_poly[6] ^ i_poly[8] ^ i_poly[9] ^ i_poly[3] ^ i_poly[5];
			shift_poly_62_10[5] = i_poly[3] ^ i_poly[1] ^ i_poly[7] ^ i_poly[9] ^ i_poly[4] ^ i_poly[6];
			shift_poly_62_10[6] = i_poly[4] ^ i_poly[2] ^ i_poly[8] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_62_10[7] = i_poly[5] ^ i_poly[3] ^ i_poly[9] ^ i_poly[0] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_62_10[8] = i_poly[6] ^ i_poly[4] ^ i_poly[1] ^ i_poly[0] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_62_10[9] = i_poly[7] ^ i_poly[5] ^ i_poly[2] ^ i_poly[1] ^ i_poly[3] ^ i_poly[8];
			shift_poly_62_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_63_10;
		input [9:0] i_poly;
		begin
			shift_poly_63_10[0] = i_poly[7] ^ i_poly[5] ^ i_poly[2] ^ i_poly[1] ^ i_poly[3] ^ i_poly[8];
			shift_poly_63_10[1] = i_poly[8] ^ i_poly[6] ^ i_poly[0] ^ i_poly[3] ^ i_poly[2] ^ i_poly[4] ^ i_poly[9];
			shift_poly_63_10[2] = i_poly[9] ^ i_poly[7] ^ i_poly[1] ^ i_poly[4] ^ i_poly[3] ^ i_poly[5];
			shift_poly_63_10[3] = i_poly[0] ^ i_poly[4] ^ i_poly[6] ^ i_poly[7] ^ i_poly[1] ^ i_poly[3];
			shift_poly_63_10[4] = i_poly[1] ^ i_poly[5] ^ i_poly[7] ^ i_poly[8] ^ i_poly[2] ^ i_poly[4];
			shift_poly_63_10[5] = i_poly[2] ^ i_poly[0] ^ i_poly[6] ^ i_poly[8] ^ i_poly[9] ^ i_poly[3] ^ i_poly[5];
			shift_poly_63_10[6] = i_poly[3] ^ i_poly[1] ^ i_poly[7] ^ i_poly[9] ^ i_poly[4] ^ i_poly[6];
			shift_poly_63_10[7] = i_poly[4] ^ i_poly[2] ^ i_poly[8] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_63_10[8] = i_poly[5] ^ i_poly[3] ^ i_poly[9] ^ i_poly[0] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_63_10[9] = i_poly[6] ^ i_poly[4] ^ i_poly[1] ^ i_poly[0] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_63_10[10] = 1'b0;
		end
	endfunction

	function automatic [10:0] shift_poly_64_10;
		input [9:0] i_poly;
		begin
			shift_poly_64_10[0] = i_poly[6] ^ i_poly[4] ^ i_poly[1] ^ i_poly[0] ^ i_poly[2] ^ i_poly[7] ^ i_poly[9];
			shift_poly_64_10[1] = i_poly[7] ^ i_poly[5] ^ i_poly[2] ^ i_poly[1] ^ i_poly[3] ^ i_poly[8];
			shift_poly_64_10[2] = i_poly[8] ^ i_poly[6] ^ i_poly[0] ^ i_poly[3] ^ i_poly[2] ^ i_poly[4] ^ i_poly[9];
			shift_poly_64_10[3] = i_poly[3] ^ i_poly[5] ^ i_poly[6] ^ i_poly[0] ^ i_poly[2];
			shift_poly_64_10[4] = i_poly[0] ^ i_poly[4] ^ i_poly[6] ^ i_poly[7] ^ i_poly[1] ^ i_poly[3];
			shift_poly_64_10[5] = i_poly[1] ^ i_poly[5] ^ i_poly[7] ^ i_poly[8] ^ i_poly[2] ^ i_poly[4];
			shift_poly_64_10[6] = i_poly[2] ^ i_poly[0] ^ i_poly[6] ^ i_poly[8] ^ i_poly[9] ^ i_poly[3] ^ i_poly[5];
			shift_poly_64_10[7] = i_poly[3] ^ i_poly[1] ^ i_poly[7] ^ i_poly[9] ^ i_poly[4] ^ i_poly[6];
			shift_poly_64_10[8] = i_poly[4] ^ i_poly[2] ^ i_poly[8] ^ i_poly[0] ^ i_poly[5] ^ i_poly[7];
			shift_poly_64_10[9] = i_poly[5] ^ i_poly[3] ^ i_poly[9] ^ i_poly[0] ^ i_poly[1] ^ i_poly[6] ^ i_poly[8];
			shift_poly_64_10[10] = 1'b0;
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
                case (code_r) // synopsys parallel_case
                    1: shift = poly_reduce_6(temp);
                    2: shift = poly_reduce_8(temp);
                    3: shift = poly_reduce_10(temp);
                    default: shift = temp;
                endcase
            end
        end
    endfunction

	function automatic [10:0] compute_d;
		input [3:0] i_mu_plus_1; // mu + 1
		input [3:0] i_lu_plus_1; // l_{mu + 1}
		input [1:0] i_index; 
		reg [3:0] cnt;
		integer i;
		begin
			cnt = i_mu_plus_1 + 1;
			compute_d = 0;
			for (i = 0; i < 5; i = i + 1) begin
				if ($signed(i[3:0]) <= $signed(i_lu_plus_1)) begin
					// S_{μ+2-i} = S_r[μ+1-i] = S_r[i_mu_plus_1 - i]
					if ($signed(i_mu_plus_1) >= $signed(i[3:0])) begin
						compute_d = compute_d ^ element_mul(delta_w[i_index][i], S_r[i_index][cnt - 1]);
						cnt = cnt - 1;
					end
				end
			end
		end
		
	endfunction

	function automatic [9:0] compare_corr;
		input [9:0] i_corr1;
		input [9:0] i_corr2;
		input [9:0] i_corr3;
		input [9:0] i_corr4;
		reg [9:0] stage1_min1;
		reg [9:0] stage1_min2;
		reg [1:0] index1;
		reg [1:0] index2;
		begin
			index1 = 2'd0;
			index2 = 2'd0;
			compare_corr = 0;
			if (i_corr1 <= i_corr2) begin
				stage1_min1 = i_corr1;
				index1 = 2'd0;
			end 
			else begin 
				stage1_min1 = i_corr2;
				index1 = 2'd1;
			end
			if (i_corr3 <= i_corr4) begin 
				stage1_min2 = i_corr3;
				index2 = 2'd2;
			end 
			else begin 
				stage1_min2 = i_corr4;
				index2 = 2'd3;
			end
			if (stage1_min1 <= stage1_min2) compare_corr = index1;
			else compare_corr = index2;
		end
		
	endfunction

endmodule

// swiss
module swiss (
    input clk,
    input rst_n,
    input [6:0] abs_data0,
    input [6:0] abs_data1,
    input [6:0] abs_data2,
    input [6:0] abs_data3,
    input [6:0] abs_data4,
    input [6:0] abs_data5,
    input [6:0] abs_data6,
    input [6:0] abs_data7,
    input [9:0] cnt,
    output [6:0] min_1,
    output [6:0] min_2,
    output [9:0] index_1,
    output [9:0] index_2
);
    integer i;
    
    reg [6:0] stage1_1_r [0:3], stage1_1_w [0:3], stage1_0_r [0:3], stage1_0_w [0:3];
    reg [9:0] index_stage1_1_r [0:3], index_stage1_1_w [0:3], index_stage1_0_r [0:3], index_stage1_0_w [0:3];
    
    reg [6:0] stage2_min_r, stage2_min_w;
    reg [9:0] index_stage2_min_r, index_stage2_min_w;
    reg [6:0] stage2_110_r, stage2_110_w;
    reg [9:0] index_stage2_110_r, index_stage2_110_w;
    reg [6:0] stage2_101_r [0:1], stage2_101_w [0:1];
    reg [9:0] index_stage2_101_r [0:1], index_stage2_101_w [0:1];
    
    reg [6:0] stage3_min1_r, stage3_min1_w, stage3_min2_r, stage3_min2_w;
    reg [9:0] index_stage3_min1_r, index_stage3_min1_w, index_stage3_min2_r, index_stage3_min2_w;

    wire [6:0] abs_data [0:7];
    assign abs_data[0] = abs_data0;
    assign abs_data[1] = abs_data1;
    assign abs_data[2] = abs_data2;
    assign abs_data[3] = abs_data3;
    assign abs_data[4] = abs_data4;
    assign abs_data[5] = abs_data5;
    assign abs_data[6] = abs_data6;
    assign abs_data[7] = abs_data7;

    always @(*) begin
        for (i = 0; i < 4; i = i + 1) begin
            if (abs_data[2*i] < abs_data[2*i+1]) begin
                stage1_1_w[i] = abs_data[2*i];
                stage1_0_w[i] = abs_data[2*i+1];
                index_stage1_1_w[i] = cnt - 2 * i;
                index_stage1_0_w[i] = cnt - 2 * i - 1;
            end else begin
                stage1_1_w[i] = abs_data[2*i+1];
                stage1_0_w[i] = abs_data[2*i];
                index_stage1_1_w[i] = cnt - 2 * i - 1;
                index_stage1_0_w[i] = cnt - 2 * i;
            end
        end
    end

    always @(*) begin
        reg [6:0] temp_stage2_11 [0:1], temp_stage2_10 [0:1], temp_stage2_01 [0:1];
        reg [9:0] temp_index_stage2_11 [0:1], temp_index_stage2_10 [0:1], temp_index_stage2_01 [0:1];

        for (i = 0; i < 2; i = i + 1) begin
            if (stage1_1_r[2*i] < stage1_1_r[2*i+1]) begin
                temp_stage2_11[i] = stage1_1_r[2*i];
                temp_stage2_10[i] = stage1_1_r[2*i+1];
                temp_index_stage2_11[i] = index_stage1_1_r[2*i];
                temp_index_stage2_10[i] = index_stage1_1_r[2*i+1];
            end else begin
                temp_stage2_11[i] = stage1_1_r[2*i+1];
                temp_stage2_10[i] = stage1_1_r[2*i];
                temp_index_stage2_11[i] = index_stage1_1_r[2*i+1];
                temp_index_stage2_10[i] = index_stage1_1_r[2*i];
            end
        end

        for (i = 0; i < 2; i = i + 1) begin
            if (stage1_0_r[2*i] < stage1_0_r[2*i+1]) begin
                temp_stage2_01[i] = stage1_0_r[2*i];
                temp_index_stage2_01[i] = index_stage1_0_r[2*i];
            end else begin
                temp_stage2_01[i] = stage1_0_r[2*i+1];
                temp_index_stage2_01[i] = index_stage1_0_r[2*i+1];
            end
        end

        if (temp_stage2_11[0] < temp_stage2_11[1]) begin
            stage2_min_w = temp_stage2_11[0];
            stage2_110_w = temp_stage2_11[1];
            index_stage2_min_w = temp_index_stage2_11[0];
            index_stage2_110_w = temp_index_stage2_11[1];
        end else begin
            stage2_min_w = temp_stage2_11[1];
            stage2_110_w = temp_stage2_11[0];
            index_stage2_min_w = temp_index_stage2_11[1];
            index_stage2_110_w = temp_index_stage2_11[0];
        end

        for (i = 0; i < 2; i = i + 1) begin
            if (temp_stage2_10[i] < temp_stage2_01[i]) begin
                stage2_101_w[i] = temp_stage2_10[i];
                index_stage2_101_w[i] = temp_index_stage2_10[i];
            end else begin
                stage2_101_w[i] = temp_stage2_01[i];
                index_stage2_101_w[i] = temp_index_stage2_01[i];
            end
        end
    end

    always @(*) begin
        reg [6:0] temp_stage3_candidate;
        reg [9:0] temp_index_stage3_candidate;

        if (stage2_101_r[0] < stage2_101_r[1]) begin
            temp_stage3_candidate = stage2_101_r[0];
            temp_index_stage3_candidate = index_stage2_101_r[0];
        end else begin
            temp_stage3_candidate = stage2_101_r[1];
            temp_index_stage3_candidate = index_stage2_101_r[1];
        end

        if (temp_stage3_candidate < stage2_110_r) begin
            stage3_min2_w = temp_stage3_candidate;
            index_stage3_min2_w = temp_index_stage3_candidate;
        end else begin
            stage3_min2_w = stage2_110_r;
            index_stage3_min2_w = index_stage2_110_r;
        end

        stage3_min1_w = stage2_min_r;
        index_stage3_min1_w = index_stage2_min_r;
    end

    always @(posedge clk) begin
        if (!rst_n) begin
            for (i = 0; i < 4; i = i + 1) begin
                stage1_1_r[i] <= 7'd0;
                stage1_0_r[i] <= 7'd0;
                index_stage1_1_r[i] <= 10'd0;
                index_stage1_0_r[i] <= 10'd0;
            end

            stage2_min_r <= 7'd0;
            index_stage2_min_r <= 10'd0;
            stage2_110_r <= 7'd0;
            index_stage2_110_r <= 10'd0;
            for (i = 0; i < 2; i = i + 1) begin
                stage2_101_r[i] <= 7'd0;
                index_stage2_101_r[i] <= 10'd0;
            end

            stage3_min1_r <= 7'd0;
            stage3_min2_r <= 7'd0;
            index_stage3_min1_r <= 10'd0;
            index_stage3_min2_r <= 10'd0;
        end
        else begin
            for (i = 0; i < 4; i = i + 1) begin
                stage1_1_r[i] <= stage1_1_w[i];
                stage1_0_r[i] <= stage1_0_w[i];
                index_stage1_1_r[i] <= index_stage1_1_w[i];
                index_stage1_0_r[i] <= index_stage1_0_w[i];
            end
            
            stage2_min_r <= stage2_min_w;
            index_stage2_min_r <= index_stage2_min_w;
            stage2_110_r <= stage2_110_w;
            index_stage2_110_r <= index_stage2_110_w;
            for (i = 0; i < 2; i = i + 1) begin
                stage2_101_r[i] <= stage2_101_w[i];
                index_stage2_101_r[i] <= index_stage2_101_w[i];
            end
            
            stage3_min1_r <= stage3_min1_w;
            stage3_min2_r <= stage3_min2_w;
            index_stage3_min1_r <= index_stage3_min1_w;
            index_stage3_min2_r <= index_stage3_min2_w;
        end
    end

    assign min_1 = stage3_min1_r;
    assign min_2 = stage3_min2_r;
    assign index_1 = index_stage3_min1_r;
    assign index_2 = index_stage3_min2_r;
endmodule