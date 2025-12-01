## 1 作业内容

①下载并安装PyTorch和Transformers
②使用Transformers库实现一个模型(任意)进行任意数量轮次的对话(最好超过两轮)，并保留历史记录;
③熟悉分词，并分别将一个句子(任何)转换为其分词形式和ID形式。
④熟悉模型的结构，并打印该模型任意层的一个(或多个)attention层的输出。

## 2 代码实现

### 2.1 多轮对话

目前下载好的PyTorch版本为2.5.1+cu121，Transformers版本4.52.3，CUDA版本 12.1 (nvcc 12.1.66)

<img src="D:\xwechat_files\wxid_uzqin1xmai4h22_a423\temp\RWTemp\2025-10\9e20f478899dc29eb19741386f9343c8\0691738425e99f9b72f6fa53c881a952.png" alt="0691738425e99f9b72f6fa53c881a952" style="zoom:50%;" />

选择用Qwen/Qwen3-0.6B模型，轻量且效果好。

```python
from transformers import AutoTokenizer, AutoModelForCausalLM
model_name="Qwen/Qwen3-0.6B"
#导入模型
tokenizer=AutoTokenizer.from_pretrained(model_name)
model=AutoModelForCausalLM.from_pretrained(model_name)

chat_history = []#全局变量记录历史记录
```

如果要实现带历史记录的输入，可以用tokenizer.apply_chat_template方法，可以将对话数据格式化为可接受的输入格式，支持Qwen的对话模板，能够自动处理历史记录。

chat_history以**列表嵌套字典**的格式存储历史记录，tokenizer.apply_chat_template可以自动将其转化为对话格式。

定义一个函数来调用模型进行对话。函数输入为User新的语句，模型输入为整个的历史记录，输出为模型的最后一次输出。

```python
from transformers import AutoTokenizer, AutoModelForCausalLM
import torch

model_name="Qwen/Qwen3-0.6B"

tokenizer=AutoTokenizer.from_pretrained(model_name)
model=AutoModelForCausalLM.from_pretrained(model_name)

chat_history=[]

def chat(user_input,max_new_tokens=256):
    global chat_history
    chat_history.append({"role": "user", "content": user_input})
    templated_input=tokenizer.apply_chat_template(
        chat_history,
        tokenize=False,
        add_generation_prompt=True 
    )
    #分词转为张量
    inputs=tokenizer(templated_input, return_tensors="pt").to(model.device)
    outputs=model.generate(
        **inputs,
        max_new_tokens=max_new_tokens,
    )
    #解码输出
    reply=tokenizer.decode(outputs[0], skip_special_tokens=True)
    #提取Assistant的新回复部分
    response=reply.split("assistant")[-1].replace(":", "").strip()
    #存入历史记录
    chat_history.append({"role": "assistant", "content": response})
    return response
```

输入对话：

```python
print("User: Hello, who are you?")
print("Assistant:", chat("Hello, who are you?"))

print("\nUser:What is the capital of China?")
print("Assistant:", chat("What is the capital of China?"))

print("\nUser:What local specialties are there here?")
print("Assistant:", chat("What local specialties are there here?"))

print("\nUser:What famous attractions are there here?")
print("Assistant:", chat("What famous attractions are there here?"))
```

输出结果如下，Qwen的输出格式带有思考过程，即\<think>，同时也测试了其他问题，偶尔会出现客观上的错误，可能是模型小的原因，没有那么大的数据量来训练。

![409baa44a7cf5c34ca67b5f8f68016ee](D:\xwechat_files\wxid_uzqin1xmai4h22_a423\temp\RWTemp\2025-10\9e20f478899dc29eb19741386f9343c8\409baa44a7cf5c34ca67b5f8f68016ee.png)

![3c7552746dfb36af0628ebe2cf7c1e7e](D:\xwechat_files\wxid_uzqin1xmai4h22_a423\temp\RWTemp\2025-10\9e20f478899dc29eb19741386f9343c8\3c7552746dfb36af0628ebe2cf7c1e7e.png)

### 2.2 分词

每个模型自带一个分词器tokenizer，可以将一句话分为各个带有id的tokens。以便进行后续的训练。

![b2df2cf0a9bd5ca5920893bd90211773](D:\xwechat_files\wxid_uzqin1xmai4h22_a423\temp\RWTemp\2025-10\9e20f478899dc29eb19741386f9343c8\b2df2cf0a9bd5ca5920893bd90211773.png)

### 2.3 Attention层输出

图为Qwen3-32B模型的结构，0.6B只是参数更小：

<img src="C:\Users\Lenovo\Downloads\75a87897aee936ee962ab22d4372494b-917df26a-1b4b-43fb-b244-2ae9a0b850e4.png" alt="75a87897aee936ee962ab22d4372494b-917df26a-1b4b-43fb-b244-2ae9a0b850e4" style="zoom: 33%;" />

据Qwen3的模型文档以及代码，Qwen3共有28层，每层共有16个注意力头。

可在outputs时设置是否输出attention层，此处输出第0层的attention层的16个attention头的输出

![9c48ad90a72f34bb095ce2eee66ece09](D:\xwechat_files\wxid_uzqin1xmai4h22_a423\temp\RWTemp\2025-10\9e20f478899dc29eb19741386f9343c8\9c48ad90a72f34bb095ce2eee66ece09.png)

<img src="D:\xwechat_files\wxid_uzqin1xmai4h22_a423\temp\RWTemp\2025-10\9e20f478899dc29eb19741386f9343c8\2c52569af49c070ecd548138b9805c7f.png" alt="2c52569af49c070ecd548138b9805c7f" style="zoom: 38%;" /><img src="D:\xwechat_files\wxid_uzqin1xmai4h22_a423\temp\RWTemp\2025-10\9e20f478899dc29eb19741386f9343c8\bba66320dba979e90f27cf3a89d14804.png" alt="bba66320dba979e90f27cf3a89d14804" style="zoom: 38%;" /><img src="D:\xwechat_files\wxid_uzqin1xmai4h22_a423\temp\RWTemp\2025-10\9e20f478899dc29eb19741386f9343c8\139d2fe4f474e9ce62fed3c7b32f975f.png" alt="139d2fe4f474e9ce62fed3c7b32f975f" style="zoom: 38%;" /><img src="D:\xwechat_files\wxid_uzqin1xmai4h22_a423\temp\RWTemp\2025-10\9e20f478899dc29eb19741386f9343c8\e9762c0687f0e85c3fb9742b17c9eb59.png" alt="e9762c0687f0e85c3fb9742b17c9eb59" style="zoom: 38%;" />