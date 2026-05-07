module async_data_mux (
    input  wire rst_n,
    input  wire clk_a,
    input  wire data_a,
    input  wire clk_b,
    input  wire data_b,
    input  wire clk_out,
    output reg  data_out,
    output reg [1:0] frame_type,
    output reg  overflow_a,
    output reg  overflow_b
);

    localparam [7:0] SYNC_A = 8'hAA;
    localparam [7:0] SYNC_B = 8'hCC;
    localparam [31:0] IDLE_FRAME = {8'b11101110, 24'h000000};

    reg [31:0] shift_a;
    reg [31:0] shift_b;
    reg [31:0] wr_data_a;
    reg [31:0] wr_data_b;
    reg        wr_en_a;
    reg        wr_en_b;

    wire [31:0] fifo_a_dout;
    wire [31:0] fifo_b_dout;
    wire        fifo_a_full;
    wire        fifo_b_full;
    wire        fifo_a_empty;
    wire        fifo_b_empty;

    reg         rd_en_a;
    reg         rd_en_b;

    reg [31:0] tx_shift;
    reg [5:0]  tx_bit_cnt;
    reg        last_sel;

    always @(posedge clk_a or negedge rst_n) begin
        if (!rst_n) begin
            shift_a  <= 32'd0;
            wr_data_a <= 32'd0;
            wr_en_a  <= 1'b0;
            overflow_a <= 1'b0;
        end else begin
            shift_a <= {shift_a[30:0], data_a};
            wr_en_a <= 1'b0;
            if ({shift_a[30:0], data_a}[31:24] == SYNC_A) begin
                if (!fifo_a_full) begin
                    wr_data_a <= {shift_a[30:0], data_a};
                    wr_en_a   <= 1'b1;
                end else begin
                    overflow_a <= 1'b1;
                end
            end
        end
    end

    always @(posedge clk_b or negedge rst_n) begin
        if (!rst_n) begin
            shift_b  <= 32'd0;
            wr_data_b <= 32'd0;
            wr_en_b  <= 1'b0;
            overflow_b <= 1'b0;
        end else begin
            shift_b <= {shift_b[30:0], data_b};
            wr_en_b <= 1'b0;
            if ({shift_b[30:0], data_b}[31:24] == SYNC_B) begin
                if (!fifo_b_full) begin
                    wr_data_b <= {shift_b[30:0], data_b};
                    wr_en_b   <= 1'b1;
                end else begin
                    overflow_b <= 1'b1;
                end
            end
        end
    end

    async_fifo #(.DATA_WIDTH(32), .ADDR_WIDTH(4)) u_fifo_a (
        .wr_clk (clk_a),
        .rd_clk (clk_out),
        .rst_n  (rst_n),
        .wr_en  (wr_en_a),
        .wr_data(wr_data_a),
        .full   (fifo_a_full),
        .rd_en  (rd_en_a),
        .rd_data(fifo_a_dout),
        .empty  (fifo_a_empty)
    );

    async_fifo #(.DATA_WIDTH(32), .ADDR_WIDTH(4)) u_fifo_b (
        .wr_clk (clk_b),
        .rd_clk (clk_out),
        .rst_n  (rst_n),
        .wr_en  (wr_en_b),
        .wr_data(wr_data_b),
        .full   (fifo_b_full),
        .rd_en  (rd_en_b),
        .rd_data(fifo_b_dout),
        .empty  (fifo_b_empty)
    );

    always @(posedge clk_out or negedge rst_n) begin
        if (!rst_n) begin
            rd_en_a   <= 1'b0;
            rd_en_b   <= 1'b0;
            tx_shift  <= {IDLE_FRAME[30:0], 1'b0};
            tx_bit_cnt <= 6'd0;
            data_out  <= 1'b0;
            frame_type <= 2'd0;
            last_sel  <= 1'b0;
        end else begin
            rd_en_a <= 1'b0;
            rd_en_b <= 1'b0;
            if (tx_bit_cnt == 6'd0) begin
                if (!fifo_a_empty && !fifo_b_empty) begin
                    if (last_sel) begin
                        rd_en_a    <= 1'b1;
                        tx_shift   <= {fifo_a_dout[30:0], 1'b0};
                        data_out   <= fifo_a_dout[31];
                        frame_type <= 2'd1;
                        last_sel   <= 1'b0;
                    end else begin
                        rd_en_b    <= 1'b1;
                        tx_shift   <= {fifo_b_dout[30:0], 1'b0};
                        data_out   <= fifo_b_dout[31];
                        frame_type <= 2'd2;
                        last_sel   <= 1'b1;
                    end
                end else if (!fifo_a_empty) begin
                    rd_en_a    <= 1'b1;
                    tx_shift   <= {fifo_a_dout[30:0], 1'b0};
                    data_out   <= fifo_a_dout[31];
                    frame_type <= 2'd1;
                    last_sel   <= 1'b0;
                end else if (!fifo_b_empty) begin
                    rd_en_b    <= 1'b1;
                    tx_shift   <= {fifo_b_dout[30:0], 1'b0};
                    data_out   <= fifo_b_dout[31];
                    frame_type <= 2'd2;
                    last_sel   <= 1'b1;
                end else begin
                    tx_shift   <= {IDLE_FRAME[30:0], 1'b0};
                    data_out   <= IDLE_FRAME[31];
                    frame_type <= 2'd0;
                end
                tx_bit_cnt <= 6'd31;
            end else begin
                data_out   <= tx_shift[31];
                tx_shift   <= {tx_shift[30:0], 1'b0};
                tx_bit_cnt <= tx_bit_cnt - 1'b1;
            end
        end
    end

endmodule
