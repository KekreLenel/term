function ImportData()

    width = 1.662557349412794

    csv_reader = CSV.File("../data/highfrequency/yieldcurve_bymaturity.csv")

    df1_real_series = csv_reader.df1_real
    df1_real_lci    = df1_real_series-width*csv_reader.df1_real_se
    df1_real_uci    = df1_real_series+width*csv_reader.df1_real_se
    df1_real        = Empirics(df1_real_series, df1_real_lci, df1_real_uci)

    df1_series = csv_reader.df1
    df1_lci    = df1_series-width*csv_reader.df1_se
    df1_uci    = df1_series+width*csv_reader.df1_se
    df1        = Empirics(df1_series, df1_lci, df1_uci)

    dy1_series = csv_reader.dy1
    dy1_lci    = dy1_series-width*csv_reader.dy1_se
    dy1_uci    = dy1_series+width*csv_reader.dy1_se
    dy1        = Empirics(dy1_series, dy1_lci, dy1_uci)

    csv_reader = CSV.File("../data/moments/tab_bymaturity.csv")

    fb_series = csv_reader.fb
    fb_lci    = fb_series-width*csv_reader.fb_se
    fb_uci    = fb_series+width*csv_reader.fb_se
    fb        = Empirics(fb_series, fb_lci, fb_uci)

    cs_series = csv_reader.csh
    cs_lci    = cs_series-width*csv_reader.csh_se
    cs_uci    = cs_series+width*csv_reader.csh_se
    cs        = Empirics(cs_series, cs_lci, cs_uci)

    csv_reader = CSV.File("../data/highfrequency/scatterplots.csv")

    year       = string.(csv_reader.year)
    year       = SubString.(year,3,4)
    month      = string.(csv_reader.month)
    yms        = month.*"/".*year
    dy1_real_hat = csv_reader.dy_real_1_hat
    df1_real_20 = csv_reader.df1_real_20
    avg_dealer_return = csv_reader.avg_dealer_return
    dy1_real_hat_hf      = csv_reader.dy_real_1_hat
    avg_dealer_return_hf = csv_reader.avg_dealer_return_hf
    yms_hf     = month.*"/".*year

    tips_issues = XLSX.readdata("../data/tipsissue/TIPS_issue.xlsx","Sheet1", "A1:BK240")

    csv_reader = CSV.File("../data/duration/duration.csv")

    tp_5yf5y   = csv_reader.dkw_realtermprem5f5[10598:13852]
    log_dur_pd = csv_reader.log_dur_pd[10598:13852]
    aggincgap  = csv_reader.aggincgap[10598:13852]
    years      = string.(csv_reader.year[10598:13852])
    year       = SubString.(years,3,4)
    month      = string.(csv_reader.month[10598:13852])
    yd         = years

    csv_reader = CSV.File("../data/qe/purchases.csv")

    qe1_10ye   = csv_reader.qe1_10ye
    year       = string.(csv_reader.year)
    year       = SubString.(year,3,4)
    month      = string.(csv_reader.month)
    ymqe       = month.*"/".*year

    df = CSV.read("../data/qe/purchasesZCP_surprise.csv",DataFrame)
    qe_purchases = Matrix(df)
    qe_purchases = qe_purchases[6:end,3:end]'

    data_empirical = RunEmpirics(fb, cs, df1_real, df1, dy1, yms, yms_hf, dy1_real_hat, dy1_real_hat_hf, df1_real_20, 
    avg_dealer_return, avg_dealer_return_hf, tips_issues, qe_purchases, tp_5yf5y, log_dur_pd, aggincgap, yd, qe1_10ye, ymqe)

    return data_empirical

end
