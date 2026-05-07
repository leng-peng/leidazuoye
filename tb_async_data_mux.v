`timescale 1ns/1ps

module tb_async_data_mux;

    reg rst_n = 1'b0;
    reg clk_a = 1'b0;
    reg clk_b = 1'b0;
    reg clk_out = 1'b0;
    reg data_a = 1'b0;
    reg data_b = 1'b0;

    wire data_out;
    wire [1:0] frame_type;
    wire overflow_a;
    wire overflow_b;

    reg [31:0] out_shift;
    integer out_bit_cnt;
    integer out_a_cnt;
    integer out_b_cnt;
    integer out_idle_cnt;
    integer sum_a_payload;
    integer sum_b_payload;

    async_data_mux dut (
        .rst_n(rst_n),
        .clk_a(clk_a),
        .data_a(data_a),
        .clk_b(clk_b),
        .data_b(data_b),
        .clk_out(clk_out),
        .data_out(data_out),
        .frame_type(frame_type),
        .overflow_a(overflow_a),
        .overflow_b(overflow_b)
    );

    always #16.666 clk_a = ~clk_a;   // 30 MHz
    always #10.000 clk_b = ~clk_b;   // 50 MHz
    always #5.682  clk_out = ~clk_out; // 88 MHz

    task send_frame_a(input [31:0] frame);
        integer i;
        begin
            for (i = 31; i >= 0; i = i - 1) begin
                @(posedge clk_a);
                data_a <= frame[i];
            end
        end
    endtask

    task send_frame_b(input [31:0] frame);
        integer i;
        begin
            for (i = 31; i >= 0; i = i - 1) begin
                @(posedge clk_b);
                data_b <= frame[i];
            end
        end
    endtask

    always @(posedge clk_out) begin
        if (!rst_n) begin
            out_shift   <= 32'd0;
            out_bit_cnt <= 0;
            out_a_cnt   <= 0;
            out_b_cnt   <= 0;
            out_idle_cnt <= 0;
            sum_a_payload <= 0;
            sum_b_payload <= 0;
        end else begin
            out_shift <= {out_shift[30:0], data_out};
            if (out_bit_cnt == 31) begin
                if ({out_shift[30:0], data_out}[31:24] == 8'hAA) begin
                    out_a_cnt <= out_a_cnt + 1;
                    sum_a_payload <= sum_a_payload + {out_shift[30:0], data_out}[23:0];
                end else if ({out_shift[30:0], data_out}[31:24] == 8'hCC) begin
                    out_b_cnt <= out_b_cnt + 1;
                    sum_b_payload <= sum_b_payload + {out_shift[30:0], data_out}[23:0];
                end else if ({out_shift[30:0], data_out}[31:24] == 8'hEE) begin
                    out_idle_cnt <= out_idle_cnt + 1;
                end
                out_bit_cnt <= 0;
            end else begin
                out_bit_cnt <= out_bit_cnt + 1;
            end
        end
    end

    initial begin
        $dumpfile("async_data_mux.vcd");
        $dumpvars(0, tb_async_data_mux);

        #100;
        rst_n = 1'b1;

        fork
            begin : stream_a
                send_frame_a(32'hAA11_2233);
                send_frame_a(32'hAA44_5566);
                send_frame_a(32'hAA77_8899);
                send_frame_a(32'hAAAB_CDEF);
                send_frame_a(32'hAA01_0203);
                send_frame_a(32'hAA04_0506);
            end
            begin : stream_b
                send_frame_b(32'hCC10_2030);
                send_frame_b(32'hCC40_5060);
                send_frame_b(32'hCC70_8090);
                send_frame_b(32'hCCA0_B0C0);
                send_frame_b(32'hCCD0_E0F0);
                send_frame_b(32'hCC12_3456);
                send_frame_b(32'hCC65_4321);
                send_frame_b(32'hCCDE_ADBE);
                send_frame_b(32'hCCEF_0011);
                send_frame_b(32'hCC22_3344);
            end
        join

        #20000;

        $display("Output frames: A=%0d, B=%0d, IDLE=%0d", out_a_cnt, out_b_cnt, out_idle_cnt);
        if (out_a_cnt != 6 || out_b_cnt != 10 || out_idle_cnt == 0) begin
            $display("TEST FAILED: frame loss detected");
            $finish_and_return(1);
        end
        if (sum_a_payload != 32'h017D_D52A || sum_b_payload != 32'h0499_DB5A) begin
            $display("TEST FAILED: payload mismatch");
            $finish_and_return(1);
        end
        if (overflow_a || overflow_b) begin
            $display("TEST FAILED: input fifo overflow");
            $finish_and_return(1);
        end

        $display("TEST PASSED");
        $finish;
    end

endmodule
