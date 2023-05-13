from comcrawl import IndexClient

client = IndexClient()

client.search("reddit.com/r/MachineLearning/*")
client.download()

first_page_html = client.results[0]["html"]