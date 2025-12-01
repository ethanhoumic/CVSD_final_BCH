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

	reg [3:0] state_r, state_w;
	reg [9:0] cnt_r, cnt_w;

	// input data
	reg [7:0] data_r [0:1023], data_w [0:1023];
	reg mode_r, mode_w;
	reg [1:0] code_r, code_w;

	// syndrome calculation
	reg [10:0] reduced_data0_r [0:1023], reduced_data0_w [0:1023];
	reg [10:0] reduced_data1_r [0:1023], reduced_data1_w [0:1023];
	reg [10:0] reduced_data2_r [0:1023], reduced_data2_w [0:1023];
	reg [10:0] reduced_data3_r [0:1023], reduced_data3_w [0:1023];
	reg [10:0] reduced_data4_r [0:1023], reduced_data4_w [0:1023];
	reg [10:0] reduced_data5_r [0:1023], reduced_data5_w [0:1023];
	reg [10:0] reduced_data6_r [0:1023], reduced_data6_w [0:1023];
	reg [10:0] reduced_data7_r [0:1023], reduced_data7_w [0:1023];
	reg [10:0] S_r [0:7], S_w [0:7];
	
	//  ber algorithm
	reg [10:0] delta_r [0:4], delta_w [0:4], delta_rho_r [0:4], delta_rho_w [0:4], temp1_w [0:4], temp2_w [0:4], temp3_w[0:4];
	reg [10:0] d_r, d_w, d_rho_r, d_rho_w; 
	reg [3:0]  l_r, l_w, l_rho_r, l_rho_w;
	reg [3:0]  rho_r, rho_w;

	integer i, j;

	// clock gating
	wire load_en     = (state_r == S_LOAD);
	wire syn_low_en  = (state_r == S_SYN && code_r != 3);
	wire syn_high_en = (state_r == S_SYN && code_r == 3);
	wire ber_en      = (state_r == S_BER || (state_r == S_SYN && state_w == S_BER));

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
		for (i = 0; i < 5; i = i + 1) begin
			delta_w[i] = delta_r[i];
			delta_rho_w[i] = delta_rho_r[i];
			temp1_w[i] = 11'b0;
			temp2_w[i] = 11'b0;
			temp3_w[i] = 11'b0;
		end
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
					code_w = code;
					mode_w = mode;
					ready_w = 1;
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
				case (code_r)
					1: begin
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
							S_w[0] = reduced_data0_w[cnt_r] ^ S_r[0];
							S_w[1] = reduced_data1_w[cnt_r] ^ S_r[1];
							S_w[2] = reduced_data2_w[cnt_r] ^ S_r[2];
							S_w[3] = reduced_data3_w[cnt_r] ^ S_r[3];
						end
						else begin
							cnt_w = 0;
							state_w = S_BER;
							delta_w[0] = 1;
							delta_rho_w[0] = 1;
							for (i = 0; i < 4; i = i + 1) $display("S%d = %b", i+1, S_r[i]);
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
						if (cnt_r < 255) begin
							cnt_w = cnt_r + 1;
							for (i = 0; i < 256; i = i + 1) begin
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
							S_w[0] = reduced_data0_w[cnt_r] ^ S_r[0];
							S_w[1] = reduced_data1_w[cnt_r] ^ S_r[1];
							S_w[2] = reduced_data2_w[cnt_r] ^ S_r[2];
							S_w[3] = reduced_data3_w[cnt_r] ^ S_r[3];
						end
						else begin
							cnt_w = 0;
							state_w = S_BER;
							delta_w[0] = 1;
							delta_rho_w[0] = 1;
							for (i = 0; i < 8; i = i + 1) $display("S%d = %b", i+1, S_r[i]);
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
						if (cnt_r < 1023) begin
							cnt_w = cnt_r + 1;
							for (i = 0; i < 1024; i = i + 1) begin
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
							S_w[0] = reduced_data0_w[cnt_r] ^ S_r[0];
							S_w[1] = reduced_data1_w[cnt_r] ^ S_r[1];
							S_w[2] = reduced_data2_w[cnt_r] ^ S_r[2];
							S_w[3] = reduced_data3_w[cnt_r] ^ S_r[3];
							S_w[4] = reduced_data4_w[cnt_r] ^ S_r[4];
							S_w[5] = reduced_data5_w[cnt_r] ^ S_r[5];
							S_w[6] = reduced_data6_w[cnt_r] ^ S_r[6];
							S_w[7] = reduced_data7_w[cnt_r] ^ S_r[7];
						end
						else begin
							cnt_w = 0;
							state_w = S_BER;
							delta_w[0] = 1;
							delta_rho_w[0] = 1;
							for (i = 0; i < 4; i = i + 1) $display("S%d = %b", i+1, S_r[i]);
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

				if ((cnt_r < 4 && code_r != 3) || (cnt_r < 8 && code_r == 3)) begin
					cnt_w = cnt_r + 1;
					if (d_r == 0) begin
						for (i = 0; i < 5; i = i + 1) delta_w[i] = delta_r[i];
						l_w = l_r;
					end
					else begin
						for (i = 0; i < 5; i = i + 1) begin
							temp1_w[i] = element_mul(d_rho_r, delta_r[i]);
						end
						for (i = 0; i < 5; i = i + 1) begin
							if ($signed(i[3:0]) - $signed(cnt_r) + $signed(rho_r) >= 0) temp2_w[i] = delta_rho_r[$signed(i[3:0]) - $signed(cnt_r) + $signed(rho_r)];
							else temp2_w[i] = 11'b0;
						end
						for (i = 0; i < 5; i = i + 1) begin
							temp3_w[i] = element_mul(d_r, temp2_w[i]);
						end
						for (i = 0; i < 5; i = i + 1) begin
							delta_w[i] = temp1_w[i] ^ temp3_w[i];
						end
						l_w = ($signed(l_r) > $signed(cnt_r) + $signed(l_rho_r) - $signed(rho_r)) ? $signed(l_r) : $signed(cnt_r) + $signed(l_rho_r) - $signed(rho_r);
					end
					if (d_r != 0 && $signed(cnt_r) - $signed(l_r) > $signed(rho_r) - $signed(l_rho_r)) begin
						rho_w = cnt_r;
						l_rho_w = l_r;
						d_rho_w = d_r;
						for (i = 0; i < 5; i = i + 1) begin
							delta_rho_w[i] = delta_r[i];
						end
					end
					d_w = compute_d(cnt_w, l_w);  // mu + 1 and l_{mu + 1}
				end
				else begin
					state_w = S_OUT;
					finish = 1;
					odata = 0;
				end
			end
			S_OUT: begin
				$display("The error correcting function has below coefficients:");
				for (i = 0; i < 5; i = i + 1) begin
					$display("%b * X^ %d", delta_r[i], i[3:0]);
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
			for (i = 0; i < 10; i = i + 1) begin
				d_r[i] <= 0;
				l_r[i] <= 0;
				delta_r[i] <= 0;
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
		end
		else begin
			cnt_r <= cnt_w;
			state_r <= state_w;
			ready_r <= ready_w;
			mode_r <= mode_w;
			code_r <= code_w;
			if (load_en) begin
				for (i = 0; i < 1024; i = i + 1) begin
					data_r[i] <= data_w[i];
				end
			end
			if (syn_low_en) begin
				for (i = 0; i < 1024; i = i + 1) begin
					reduced_data0_r[i] <= reduced_data0_w[i];
					reduced_data1_r[i] <= reduced_data1_w[i];
					reduced_data2_r[i] <= reduced_data2_w[i];
					reduced_data3_r[i] <= reduced_data3_w[i];
				end
				for (i = 0; i < 4; i = i + 1) begin
					S_r[i] <= S_w[i];
				end
			end
			if (syn_high_en) begin
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
		end
	end

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

	function automatic [10:0] element_mul;
		input [10:0] i_element1;
		input [10:0] i_element2;
		reg [10:0] shift;
		reg [10:0] temp;
		reg [3:0] deg;
		reg [3:0] cnt;
		integer i;
		begin
			shift = i_element1;
			element_mul = 0;
			cnt = 0;
			temp = 0;
			deg = get_degree(i_element2);
			for (i = 0; i < 11; i = i + 1) begin
				if (cnt <= deg) begin
					if (i_element2[i]) element_mul = element_mul ^ shift;
					temp = shift << 1;
					case (code_r)
						1: shift = poly_reduce_6(temp);
						2: shift = poly_reduce_8(temp);
						3: shift = poly_reduce_10(temp);
						default: shift = temp;
					endcase
					cnt = cnt + 1;
				end
				else begin
					
				end
			end
		end
	endfunction
    
    function automatic [3:0] get_degree;
        input [10:0] i_poly;
        integer i;
        reg stop;
        begin
            get_degree = 4'd10;
            stop = 0;
            for (i = 10; i >= 0; i = i - 1) begin
                if (!i_poly[i] && !stop) get_degree = (get_degree == 0) ? 0 : get_degree - 1;
                else if (i_poly[i] && !stop) stop = 1;
                else begin
                    
                end
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
