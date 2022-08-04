sed '1d' temp.detail.xls |awk -F "\t" '{print $3}'|sort|uniq|perl code_get_url.pl
