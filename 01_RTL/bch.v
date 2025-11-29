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

	reg finish;
	reg [9:0] odata;
	reg ready_r, ready_w;
	reg mode_r, mode_w;

	localparam S_IDLE = 0;
	localparam S_LOAD = 1;
	localparam S_SYN  = 2;
	localparam S_BER  = 3;
	localparam S_CHI  = 4;
	localparam S_CORR = 5;
	localparam S_OUT  = 6;

	reg [7:0] data_r [0:1023], data_w [0:1023];
	reg [10:0] reduced_data0_r [0:1023], reduced_data0_w [0:1023];
	reg [10:0] reduced_data1_r [0:1023], reduced_data1_w [0:1023];
	reg [10:0] reduced_data2_r [0:1023], reduced_data2_w [0:1023];
	reg [10:0] reduced_data3_r [0:1023], reduced_data3_w [0:1023];
	reg [10:0] reduced_data4_r [0:1023], reduced_data4_w [0:1023];
	reg [10:0] reduced_data5_r [0:1023], reduced_data5_w [0:1023];
	reg [10:0] reduced_data6_r [0:1023], reduced_data6_w [0:1023];
	reg [10:0] reduced_data7_r [0:1023], reduced_data7_w [0:1023];
	reg [10:0] S_r [0:7], S_w [0:7];
	reg [3:0] state_r, state_w;
	reg [9:0] cnt_r, cnt_w;

	reg [9:0] cnt_max_w, cnt_max_r;
	reg [2:0] t_w;

	integer i, j;

	// clock gating
	wire load_en = (state_r == S_LOAD);
	wire syn_en  = (state_r == S_SYN);

	assign ready = ready_r;

	always @(*) begin
		state_w = state_r;
		cnt_w = cnt_r;
		ready_w = ready_r;
		cnt_max_w = cnt_max_r;
		for (i = 0; i < 1024; i = i + 1) begin
			data_w[i] = data_r[i];
			reduced_data0_w[i] = reduced_data0_r[i];
			reduced_data1_w[i] = reduced_data1_r[i];
			reduced_data2_w[i] = reduced_data2_r[i];
			reduced_data3_w[i] = reduced_data3_r[i];
			reduced_data4_w[i] = reduced_data4_r[i];
			reduced_data5_w[i] = reduced_data5_r[i];
			reduced_data6_w[i] = reduced_data6_r[i];
			reduced_data7_w[i] = reduced_data7_r[i];
		end 
		for (i = 0; i < 8; i = i + 1) S_w[i] = S_r[i];
		case (state_r)
			S_IDLE: begin
				if (set && !ready_r) begin
					case (code)
						1: cnt_w = 63;
						2: cnt_w = 255;
						3: cnt_w = 1023;
						default: cnt_w = 1023;
					endcase
					mode_w = mode;
					ready_w = 1;
					cnt_max_w = cnt_w;
					state_w = S_LOAD;
				end
				// else if (ready_r) state_w = S_LOAD;
				else state_w = S_IDLE;
			end  
			S_LOAD: begin
				if (!mode_r) begin // hard decision
					data_w[cnt_r - 7] = ($signed(idata[7:0])   >= 0) ? 0 : 1;
					data_w[cnt_r - 6] = ($signed(idata[15:8])  >= 0) ? 0 : 1;
					data_w[cnt_r - 5] = ($signed(idata[23:16]) >= 0) ? 0 : 1;
					data_w[cnt_r - 4] = ($signed(idata[31:24]) >= 0) ? 0 : 1;
					data_w[cnt_r - 3] = ($signed(idata[39:32]) >= 0) ? 0 : 1;
					data_w[cnt_r - 2] = ($signed(idata[47:40]) >= 0) ? 0 : 1;
					data_w[cnt_r - 1] = ($signed(idata[55:48]) >= 0) ? 0 : 1;
					data_w[cnt_r]     = ($signed(idata[63:56]) >= 0) ? 0 : 1;
				end
				else begin // soft decision
					data_w[cnt_r - 7] = idata[7:0];
					data_w[cnt_r - 6] = idata[15:8];
					data_w[cnt_r - 5] = idata[23:16];
					data_w[cnt_r - 4] = idata[31:24];
					data_w[cnt_r - 3] = idata[39:32];
					data_w[cnt_r - 2] = idata[47:40];
					data_w[cnt_r - 1] = idata[55:48];
					data_w[cnt_r]     = idata[63:56];
				end
				if (cnt_r == 7) begin
					ready_w = 0;
					state_w = S_SYN;
					cnt_w = 0;
				end
				else begin
					cnt_w = cnt_r - 8;
				end
			end
			S_SYN: begin
				case (cnt_max_r)
					63: begin
						if (cnt_r < 63) begin
							cnt_w = cnt_r + 1;
							for (i = 0; i < 64; i = i + 1) begin
								if (cnt_r == 0) begin
									reduced_data0_w[i] = {3'b0, data_r[i]};
									reduced_data1_w[i] = {3'b0, data_r[i]};
									reduced_data2_w[i] = {3'b0, data_r[i]};
									reduced_data3_w[i] = {3'b0, data_r[i]};
								end
								else begin
									for (j = 0; j < 1; j = j + 1) reduced_data0_w[i] = poly_reduce_6({reduced_data0_w[i][9:0], 1'b0});
									for (j = 0; j < 2; j = j + 1) reduced_data1_w[i] = poly_reduce_6({reduced_data1_w[i][9:0], 1'b0});
									for (j = 0; j < 3; j = j + 1) reduced_data2_w[i] = poly_reduce_6({reduced_data2_w[i][9:0], 1'b0});
									for (j = 0; j < 4; j = j + 1) reduced_data3_w[i] = poly_reduce_6({reduced_data3_w[i][9:0], 1'b0});
								end
							end
							S_w[0] = poly_sum(reduced_data0_w[cnt_r], S_r[0]);
							S_w[1] = poly_sum(reduced_data1_w[cnt_r], S_r[1]);
							S_w[2] = poly_sum(reduced_data2_w[cnt_r], S_r[2]);
							S_w[3] = poly_sum(reduced_data3_w[cnt_r], S_r[3]);
						end
						else begin
							cnt_w = 0;
							// state_w = S_BER;
							state_w = S_OUT;
							finish = 1;
							odata = 0;
						end
					end 
					255: begin
						if (cnt_r < 255) begin
							cnt_w = cnt_r + 1;
							for (i = 0; i < 64; i = i + 1) begin
								if (cnt_r == 0) begin
									reduced_data0_w[i] = {3'b0, data_r[i]};
									reduced_data1_w[i] = {3'b0, data_r[i]};
									reduced_data2_w[i] = {3'b0, data_r[i]};
									reduced_data3_w[i] = {3'b0, data_r[i]};
								end
								else begin
									for (j = 0; j < 1; j = j + 1) reduced_data0_w[i] = poly_reduce_8({reduced_data0_w[i][9:0], 1'b0});
									for (j = 0; j < 2; j = j + 1) reduced_data1_w[i] = poly_reduce_8({reduced_data1_w[i][9:0], 1'b0});
									for (j = 0; j < 3; j = j + 1) reduced_data2_w[i] = poly_reduce_8({reduced_data2_w[i][9:0], 1'b0});
									for (j = 0; j < 4; j = j + 1) reduced_data3_w[i] = poly_reduce_8({reduced_data3_w[i][9:0], 1'b0});
								end
							end
							S_w[0] = poly_sum(reduced_data0_w[cnt_r], S_r[0]);
							S_w[1] = poly_sum(reduced_data1_w[cnt_r], S_r[1]);
							S_w[2] = poly_sum(reduced_data2_w[cnt_r], S_r[2]);
							S_w[3] = poly_sum(reduced_data3_w[cnt_r], S_r[3]);
						end
						else begin
							cnt_w = 0;
							// state_w = S_BER;
							state_w = S_OUT;
							finish = 1;
							odata = 0;
						end
					end
					1023: begin
						if (cnt_r < 1023) begin
							cnt_w = cnt_r + 1;
							for (i = 0; i < 64; i = i + 1) begin
								if (cnt_r == 0) begin
									reduced_data0_w[i] = {3'b0, data_r[i]};
									reduced_data1_w[i] = {3'b0, data_r[i]};
									reduced_data2_w[i] = {3'b0, data_r[i]};
									reduced_data3_w[i] = {3'b0, data_r[i]};
									reduced_data4_w[i] = {3'b0, data_r[i]};
									reduced_data5_w[i] = {3'b0, data_r[i]};
									reduced_data6_w[i] = {3'b0, data_r[i]};
									reduced_data7_w[i] = {3'b0, data_r[i]};
								end
								else begin
									for (j = 0; j < 1; j = j + 1) reduced_data0_w[i] = poly_reduce_10({reduced_data0_w[i][9:0], 1'b0});
									for (j = 0; j < 2; j = j + 1) reduced_data1_w[i] = poly_reduce_10({reduced_data1_w[i][9:0], 1'b0});
									for (j = 0; j < 3; j = j + 1) reduced_data2_w[i] = poly_reduce_10({reduced_data2_w[i][9:0], 1'b0});
									for (j = 0; j < 4; j = j + 1) reduced_data3_w[i] = poly_reduce_10({reduced_data3_w[i][9:0], 1'b0});
									for (j = 0; j < 5; j = j + 1) reduced_data4_w[i] = poly_reduce_10({reduced_data4_w[i][9:0], 1'b0});
									for (j = 0; j < 6; j = j + 1) reduced_data5_w[i] = poly_reduce_10({reduced_data5_w[i][9:0], 1'b0});
									for (j = 0; j < 7; j = j + 1) reduced_data6_w[i] = poly_reduce_10({reduced_data6_w[i][9:0], 1'b0});
									for (j = 0; j < 8; j = j + 1) reduced_data7_w[i] = poly_reduce_10({reduced_data7_w[i][9:0], 1'b0});
								end
							end
							S_w[0] = poly_sum(reduced_data0_w[cnt_r], S_r[0]);
							S_w[1] = poly_sum(reduced_data1_w[cnt_r], S_r[1]);
							S_w[2] = poly_sum(reduced_data2_w[cnt_r], S_r[2]);
							S_w[3] = poly_sum(reduced_data3_w[cnt_r], S_r[3]);
							S_w[4] = poly_sum(reduced_data4_w[cnt_r], S_r[4]);
							S_w[5] = poly_sum(reduced_data5_w[cnt_r], S_r[5]);
							S_w[6] = poly_sum(reduced_data6_w[cnt_r], S_r[6]);
							S_w[7] = poly_sum(reduced_data7_w[cnt_r], S_r[7]);
						end
						else begin
							cnt_w = 0;
							// state_w = S_BER;
							state_w = S_OUT;
							finish = 1;
							odata = 0;
						end
					end
				endcase
			end
			S_BER: begin
				
			end
			S_OUT: begin
				for (i = 0; i < 8; i = i + 1) begin
					$display("S_%d = %b", i+1, S_r[i]);
				end
				state_w = S_IDLE;
				if (cnt_r == 6) begin
					finish = 0;
					state_w = S_IDLE;
				end
				cnt_w = cnt_r + 1;
			end
		endcase
	end

	always @(posedge clk or negedge rstn) begin
		if (!rstn) begin
			state_r <= S_IDLE;
			for (i = 0; i < 1024; i = i + 1) begin
				data_r[i] <= 0;
				reduced_data0_r[i] <= 0;
				reduced_data1_r[i] <= 0;
				reduced_data2_r[i] <= 0;
				reduced_data3_r[i] <= 0;
				reduced_data4_r[i] <= 0;
				reduced_data5_r[i] <= 0;
				reduced_data6_r[i] <= 0;
				reduced_data7_r[i] <= 0;
			end
			for (i = 0; i < 8; i = i + 1) begin
				S_r[i] <= 0;
			end
			cnt_r <= 0;
			ready_r <= 0;
			cnt_max_r <= 0;
			mode_r <= 0;
		end
		else begin
			cnt_r <= cnt_w;
			state_r <= state_w;
			ready_r <= ready_w;
			cnt_max_r <= cnt_max_w;
			mode_r <= mode_w;
			if (load_en) begin
				for (i = 0; i < 1024; i = i + 1) begin
					data_r[i] <= data_w[i];
				end
			end
			if (syn_en) begin
				for (i = 0; i < 1024; i = i + 1) begin
					reduced_data0_r[i] <= reduced_data0_w[i];
					reduced_data1_r[i] <= reduced_data1_w[i];
					reduced_data2_r[i] <= reduced_data2_w[i];
					reduced_data3_r[i] <= reduced_data3_w[i];
					reduced_data4_r[i] <= reduced_data4_w[i];
					reduced_data5_r[i] <= reduced_data5_w[i];
					reduced_data6_r[i] <= reduced_data6_w[i];
					reduced_data7_r[i] <= reduced_data7_w[i];
				end
				for (i = 0; i < 8; i = i + 1) begin
					S_r[i] <= S_w[i];
				end
			end
		end
	end

	function automatic [10:0] poly_sum;
		input [10:0] i_poly_1;
		input [10:0] i_poly_2;
		integer i;
		begin
			poly_sum = 0;
			for (i = 0; i < 11; i = i + 1) begin
				poly_sum[i] = i_poly_1[i] ^ i_poly_2[i];
			end
		end
		
	endfunction

	function automatic [10:0] poly_reduce_6;
		input [10:0] i_poly;
		begin
			poly_reduce_6 = i_poly;
			if (i_poly[6]) begin
				poly_reduce_6[6] = 0;
				poly_reduce_6[0] = !poly_reduce_6[0];
				poly_reduce_6[1] = !poly_reduce_6[1];
			end
		end
		
	endfunction

	function automatic [10:0] poly_reduce_8;
		input [10:0] i_poly;
		begin
			poly_reduce_8 = i_poly;
			if (i_poly[8]) begin
				poly_reduce_8[8] = 0;
				poly_reduce_8[0] = !poly_reduce_8[0];
				poly_reduce_8[2] = !poly_reduce_8[2];
				poly_reduce_8[3] = !poly_reduce_8[3];
				poly_reduce_8[4] = !poly_reduce_8[4];
			end
		end
		
	endfunction

	function automatic [10:0] poly_reduce_10;
		input [10:0] i_poly;
		begin
			poly_reduce_10 = i_poly;
			if (i_poly[10]) begin
				poly_reduce_10[10] = 0;
				poly_reduce_10[0] = !poly_reduce_10[0];
				poly_reduce_10[3] = !poly_reduce_10[3];
			end
		end
		
	endfunction
	

endmodule
