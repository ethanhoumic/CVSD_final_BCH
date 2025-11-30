`timescale 1ns/1ps

module extend_euclidean (
    input         i_clk,
    input         i_rst_n,
    input         i_start,
    input [1:0]   i_code,
    input [10:0]  i_poly,
    output        o_done,
    output [10:0] o_poly
);
    
    // state
    localparam S_IDLE = 0;
    localparam S_DIV  = 1;
    localparam S_DONE = 2;

    // primitive poly
    localparam P_6  = 11'b000_0100_0011;
    localparam P_8  = 11'b001_0001_1101;
    localparam P_10 = 11'b100_0000_1001;

    reg state_r, state_w;
    reg [10:0] r0_r, r1_r, r0_w, r1_w, t0_r, t1_r, t0_w, t1_w;
    reg [1:0] code_r, code_w;

    always @(*) begin
        state_w = state_r;
        code_w = code_r;
        case (state_r)
            S_IDLE: begin
                if (i_start) begin
                    state_w = S_DIV;
                    case (i_code)
                        0:  r0_w = P_6;
                        1:  r0_w = P_8;
                        2:  r0_w = P_10;
                        default: r0_w = P_6;
                    endcase
                    r1_w = i_poly;
                    t0_w = 0;
                    t1_w = 1;
                end
            end 
            S_DIV: begin
                
            end
            default: begin
                
            end
        endcase
    end

    always @(posedge i_clk or negedge i_rst_n) begin
        if (!i_rst_n) begin
            state_r <= S_IDLE;
            r0_r <= 0;
            r1_r <= 0;
            t0_r <= 0;
            t1_r <= 0;
        end
        else begin
            state_r <= state_w;
            r0_r <= r0_w;
            r1_r <= r1_w;
            t0_r <= t0_w;
            t1_r <= t1_w;
        end
    end

endmodule

module gf_mul_unit (
    input         i_clk,
    input         i_rst_n,
    input         i_start,
    input  [1:0]  i_code,
    input  [10:0] i_poly1,
    input  [10:0] i_poly2,
    output        o_done,
    output [10:0] o_poly
);

    localparam S_IDLE = 0;
    localparam S_CAL  = 1;

    reg [10:0] poly1_r, poly2_r, poly1_w, poly2_w, shift_r, shift_w;
    reg [10:0] result_r, result_w;
    reg [3:0]  cnt_r, cnt_w;
    reg [1:0]  code_r, code_w;
    reg [10:0] temp_w;
    reg [3:0]  deg_r, deg_w;
    reg done_r, done_w;
    reg state_r, state_w;

    assign o_done = done_r;
    assign o_poly = result_r;

    always @(*) begin
        poly1_w = poly1_r;
        poly2_w = poly2_r;
        done_w=  done_r;
        state_w = state_r;
        result_w = result_r;
        cnt_w = cnt_r;
        shift_w = shift_r;
        code_w = code_r;
        temp_w = 0;
        deg_w = deg_r;
        case (state_r)
            S_IDLE: begin
                done_w = 0;
                result_w = 0;
                cnt_w = 0;
                if (i_start) begin
                    poly1_w = i_poly1;
                    poly2_w = i_poly2;
                    state_w = S_CAL;
                    shift_w = i_poly1;
                    code_w  = i_code;
                    deg_w   = get_degree(i_poly2);
                end
            end 
            S_CAL: begin
                cnt_w = cnt_r + 1;
                if (poly2_r[cnt_r]) begin
                    result_w = result_r ^ shift_r;
                end
                temp_w = shift_r << 1;
                case (code_r) // synopsys parallel case
                    0: begin
                        if (temp_w[6]) shift_w = {5'b0, temp_w[5:2], ~temp_w[1:0]};
                        else shift_w = temp_w;
                    end
                    1: begin
                        if (temp_w[8]) shift_w = {3'b0, temp_w[7:5], ~temp_w[4:2], temp_w[1], !temp_w[0]};
                        else shift_w = temp_w;
                    end
                    2: begin
                        if (temp_w[10]) shift_w = {1'b0, temp_w[9:4], !temp_w[3], temp_w[2:1], !temp_w[0]};
                        else shift_w = temp_w;
                    end
                    default: begin
                        shift_w = temp_w;
                    end
                endcase
                if (cnt_r == deg_r) begin
                    state_w = S_IDLE;
                    done_w = 1;
                    shift_w = 0;
                    cnt_w = 0;
                end
            end
        endcase
    end

    always @(posedge i_clk or negedge i_rst_n) begin
        if (!i_rst_n) begin
            poly1_r  <= 0;
            poly2_r  <= 0;
            done_r   <= 0;
            state_r  <= S_IDLE;
            result_r <= 0;
            cnt_r    <= 0;
            shift_r  <= 0;
            deg_r    <= 0;
            code_r   <= 0;
        end
        else begin
            state_r  <= state_w;
            poly1_r  <= poly1_w;
            poly2_r  <= poly2_w;
            done_r   <= done_w;
            result_r <= result_w;
            cnt_r    <= cnt_w;
            shift_r  <= shift_w;
            deg_r    <= deg_w;
            code_r   <= code_w;
        end
    end

    function automatic [3:0] get_degree;
        input [10:0] i_poly;
        integer i;
        reg stop;
        begin
            get_degree = 4'd10;
            stop = 0;
            for (i = 10; i >= 0; i = i - 1) begin
                if (!i_poly[i] && !stop) get_degree = get_degree - 1;
                else if (i_poly[i] && !stop) stop = 1;
                else begin
                    
                end
            end
        end
        
    endfunction

endmodule