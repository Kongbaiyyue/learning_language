# imports
import ast  # for converting embeddings saved as strings back to arrays
import openai  # for calling the OpenAI API
import pandas as pd  # for storing text and embeddings data
from scipy import spatial  # for calculating vector similarities for search
import os


# models
EMBEDDING_MODEL = "text-embedding-ada-002"
GPT_MODEL = "gpt-3.5-turbo"
openai.api_key = "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"

query = 'Which athletes won the gold medal in curling at the 2022 Winter Olympics?'

response = openai.ChatCompletion.create(
    messages=[
        {'role': 'system', 'content': 'You answer questions about the 2022 Winter Olympics.'},
        {'role': 'user', 'content': query},
    ],
    model=GPT_MODEL,
    temperature=0,
)

print(response['choices'][0]['message']['content'])