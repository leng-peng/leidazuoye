module async_fifo #(
    parameter DATA_WIDTH = 32,
    parameter ADDR_WIDTH = 4
) (
    input  wire                  wr_clk,
    input  wire                  rd_clk,
    input  wire                  rst_n,
    input  wire                  wr_en,
    input  wire [DATA_WIDTH-1:0] wr_data,
    output wire                  full,
    input  wire                  rd_en,
    output wire [DATA_WIDTH-1:0] rd_data,
    output wire                  empty
);

    localparam DEPTH = (1 << ADDR_WIDTH);

    reg [DATA_WIDTH-1:0] mem [0:DEPTH-1];

    reg [ADDR_WIDTH:0] wptr_bin;
    reg [ADDR_WIDTH:0] wptr_gray;
    reg [ADDR_WIDTH:0] rptr_bin;
    reg [ADDR_WIDTH:0] rptr_gray;

    reg [ADDR_WIDTH:0] wptr_gray_sync1, wptr_gray_sync2;
    reg [ADDR_WIDTH:0] rptr_gray_sync1, rptr_gray_sync2;

    wire [ADDR_WIDTH:0] wptr_bin_next;
    wire [ADDR_WIDTH:0] wptr_gray_next;
    wire [ADDR_WIDTH:0] rptr_bin_next;
    wire [ADDR_WIDTH:0] rptr_gray_next;

    assign wptr_bin_next  = wptr_bin + ((wr_en && !full) ? 1'b1 : 1'b0);
    assign wptr_gray_next = (wptr_bin_next >> 1) ^ wptr_bin_next;
    assign rptr_bin_next  = rptr_bin + ((rd_en && !empty) ? 1'b1 : 1'b0);
    assign rptr_gray_next = (rptr_bin_next >> 1) ^ rptr_bin_next;

    assign full  = (wptr_gray_next == {~rptr_gray_sync2[ADDR_WIDTH:ADDR_WIDTH-1], rptr_gray_sync2[ADDR_WIDTH-2:0]});
    assign empty = (rptr_gray == wptr_gray_sync2);

    assign rd_data = mem[rptr_bin[ADDR_WIDTH-1:0]];

    always @(posedge wr_clk or negedge rst_n) begin
        if (!rst_n) begin
            wptr_bin  <= {ADDR_WIDTH+1{1'b0}};
            wptr_gray <= {ADDR_WIDTH+1{1'b0}};
        end else begin
            if (wr_en && !full) begin
                mem[wptr_bin[ADDR_WIDTH-1:0]] <= wr_data;
                wptr_bin  <= wptr_bin_next;
                wptr_gray <= wptr_gray_next;
            end
        end
    end

    always @(posedge rd_clk or negedge rst_n) begin
        if (!rst_n) begin
            rptr_bin  <= {ADDR_WIDTH+1{1'b0}};
            rptr_gray <= {ADDR_WIDTH+1{1'b0}};
        end else begin
            if (rd_en && !empty) begin
                rptr_bin  <= rptr_bin_next;
                rptr_gray <= rptr_gray_next;
            end
        end
    end

    always @(posedge wr_clk or negedge rst_n) begin
        if (!rst_n) begin
            rptr_gray_sync1 <= {ADDR_WIDTH+1{1'b0}};
            rptr_gray_sync2 <= {ADDR_WIDTH+1{1'b0}};
        end else begin
            rptr_gray_sync1 <= rptr_gray;
            rptr_gray_sync2 <= rptr_gray_sync1;
        end
    end

    always @(posedge rd_clk or negedge rst_n) begin
        if (!rst_n) begin
            wptr_gray_sync1 <= {ADDR_WIDTH+1{1'b0}};
            wptr_gray_sync2 <= {ADDR_WIDTH+1{1'b0}};
        end else begin
            wptr_gray_sync1 <= wptr_gray;
            wptr_gray_sync2 <= wptr_gray_sync1;
        end
    end

endmodule
