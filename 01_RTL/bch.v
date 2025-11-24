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

	localparam S_IDLE = 0;
	localparam S_LOAD = 1;
	localparam S_SYN  = 2;
	localparam S_BER  = 3;
	localparam S_CHI  = 4;
	localparam S_CORR = 5;
	localparam S_OUT  = 6;

	reg [7:0] data_r [0:1023], data_w [0:1023];
	reg [10:0] reduced_data_r [0:1023], reduced_data_w [0:1023];
	reg [10:0] S_r [0:7], S_w [0:7];
	reg [3:0] state_r, state_w;
	reg [9:0] cnt_r, cnt_w;

	reg [9:0] cnt_max_w, cnt_max_r;
	reg [2:0] t_w;

	integer i;

	// clock gating
	wire load_en = (state_r == S_LOAD);
	wire syn_en  = (state_r == S_SYN);

	assign ready = ready_r;

	always @(*) begin
		state_w = state_r;
		cnt_w = cnt_r;
		ready_w = ready_r;
		cnt_max_w = cnt_max_r;
		for (i = 0; i < 1024; i = i + 1) data_w[i] = data_r[i];
		for (i = 0; i < 1024; i = i + 1) reduced_data_w[i] = reduced_data_r[i];
		case (state_r)
			S_IDLE: begin
				if (set && !ready_r) begin
					case (code)
						1: cnt_w = 63;
						2: cnt_w = 255;
						3: cnt_w = 1023;
						default: cnt_w = 1023;
					endcase
					ready_w = 1;
					cnt_max_w = cnt_w;
					state_w = S_LOAD;
				end
				// else if (ready_r) state_w = S_LOAD;
				else state_w = S_IDLE;
			end  
			S_LOAD: begin
				if (!mode) begin // hard decision
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
									reduced_data_w[i] = {3'b0, data_r[i]};
								end
								else begin
									reduced_data_w[i] = (i[9:0] >= cnt_r) ? {reduced_data_r[i][9:0], 1'b0} : reduced_data_r[i];
								end
								if (reduced_data_w[i][6]) begin
									// $display("reduce data_r[%d] at cnt = %d", i, cnt_r);
									reduced_data_w[i][6] = 0;
									reduced_data_w[i][0] = !reduced_data_w[i][0]; // xor 1
									reduced_data_w[i][1] = !reduced_data_w[i][1]; // xor 1
								end
							end
							for (i = 0; i < 63; i = i + 1) begin
								data_w[i] = data_r[i+1];
							end
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
							for (i = 0; i < 256; i = i + 1) begin
								if (cnt_r == 0) begin
									reduced_data_w[i] = {3'b0, data_r[i]};
								end
								else begin
									reduced_data_w[i] = (i[9:0] >= cnt_r) ? {reduced_data_r[i][9:0], 1'b0} : reduced_data_r[i];
								end
								if (reduced_data_w[i][8]) begin
									reduced_data_w[i][8] = 0;
									reduced_data_w[i][0] = !reduced_data_w[i][0]; // xor 1
									reduced_data_w[i][2] = !reduced_data_w[i][2]; // xor 1
									reduced_data_w[i][3] = !reduced_data_w[i][3]; // xor 1
									reduced_data_w[i][4] = !reduced_data_w[i][4]; // xor 1
								end
							end
							for (i = 0; i < 255; i = i + 1) begin
								data_w[i] = data_r[i+1];
							end
							cnt_w = cnt_r + 1;
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
							for (i = 0; i < 1024; i = i + 1) begin
								if (cnt_r == 0) begin
									reduced_data_w[i] = {3'b0, data_r[i]};
								end
								else begin
									reduced_data_w[i] = (i[9:0] >= cnt_r) ? {reduced_data_r[i][9:0], 1'b0} : reduced_data_r[i];
								end
								if (reduced_data_w[i][10]) begin
									reduced_data_w[i][10] = 0;
									reduced_data_w[i][0] = !reduced_data_w[i][0]; // xor 1
									reduced_data_w[i][3] = !reduced_data_w[i][3]; // xor 1
								end
							end
							for (i = 0; i < 1023; i = i + 1) begin
								data_w[i] = data_r[i+1];
							end
							cnt_w = cnt_r + 1;
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
				for (i = 0; i < 64; i = i + 1) begin
					// $display("reduced_data_r[%d] = %b", i, reduced_data_r[i]);
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
			end
			for (i = 0; i < 1024; i = i + 1) begin
				reduced_data_r[i] <= 0;
			end
			cnt_r <= 0;
			ready_r <= 0;
			cnt_max_r <= 0;
		end
		else begin
			cnt_r <= cnt_w;
			state_r <= state_w;
			ready_r <= ready_w;
			cnt_max_r <= cnt_max_w;
			if (load_en) begin
				for (i = 0; i < 1024; i = i + 1) begin
					data_r[i] <= data_w[i];
				end
			end
			if (syn_en) begin
				for (i = 0; i < 1024; i = i + 1) begin
					reduced_data_r[i] <= reduced_data_w[i];
				end
			end
		end
	end

endmodule
