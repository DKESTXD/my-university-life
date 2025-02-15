// components/detail_page/detail_page.js
const app=getApp();
Component({

  /**
   * 组件的属性列表
   */
  properties: {
    getthing:{
      type:Object,
      value:{}
    }
  },

  /**
   * 组件的初始数据
   */
  data: {
    "showModalViewCancel":false,
    "titleValue":"",
    "textValue":"",
    "inputLength":0,
    "star1":true,"star2":true,"star3":true,"star4":true,"star5":true,
    "urgent_statue":"",
    "urgent_level":1,
    "reminder":true,
    "settingValue":"仅仅本次",
    "timeInput":"2025-1-30-00-00",
    "value":[0,0,29,"00","00"],
    "hour":"00",
    "minute":"00",
    "year":2025,
    "month":1,
    "day":30,
    "_id":"",
    "setTime":"",
    "reminderold":true,
    "titleValueold":"",
    "urgent_levelold":1,
    "textValueold":"",
    "timeInputold":"",
    "settingValueold":"",
  },
  attached(){
    const level=this.properties.getthing.urgent_level;
    switch (level){
      case 1:
        this.setData({
          star1:true,star2:false,star3:false,star4:false,star5:false,urgent_statue:'以后的事以后再说',
          urgent_level:1
        });
        break;
      case 2:
        this.setData({
          star1:true,star2:true,star3:false,star4:false,star5:false,urgent_statue:'轻轻松 时间还很充足',
          urgent_level:2
        });
        break;
      case 3:
        this.setData({
          star1:true,star2:true,star3:true,star4:false,star5:false,urgent_statue:'你意识到你不能再拖了',
          urgent_level:3
        });
        break;
      case 4:
        this.setData({
          star1:true,star2:true,star3:true,star4:true,star5:false,urgent_statue:'非常重要 但是似乎会有比它更重要的事',urgent_level:4
        });
        break;
      case 5:
        this.setData({
          star1:true,star2:true,star3:true,star4:true,star5:true,urgent_statue:'🔥🔥🔥十万火急🔥🔥🔥',
          urgent_level:5
        });
        break;
    };
    if(this.properties.getthing.timeInput){
      const res=this.properties.getthing.timeInput.split("-")
      const year = res[0].substring(2)
      const newRes = [ year-25, res[1]-1, res[2]-1,res[3],res[4]]
      this.setData({
        value:newRes,
        year:Number(res[0]),
        month:Number(res[1]),
        day:Number(res[2]),
        hour:res[3],
        minute:res[4]
      })
    }
    this.setData({
      titleValue:this.properties.getthing.titleValue,
      titleValueold:this.properties.getthing.titleValue,
      textValue:this.properties.getthing.textValue,
      textValueold:this.properties.getthing.textValue,
      inputLength:this.properties.getthing.textValue.length,
      urgent_level:this.properties.getthing.urgent_level,
      urgent_levelold:this.properties.getthing.urgent_level,
      reminder:this.properties.getthing.reminder,
      reminderold:this.properties.getthing.reminder,
      settingValue:this.properties.getthing.settingValue,
      settingValueold:this.properties.getthing.settingValue,
      timeInput:this.properties.getthing.timeInput,
      timeInputold:this.properties.getthing.timeInput,
      _id:this.properties.getthing._id,
      setTime:this.properties.getthing.setTime,
    })
  },
  /**
   * 组件的方法列表
   */
  methods: {
    testto:function(){
      wx.navigateTo({
        url: `../testpage/testpage`
      })
    },
    edit:function() {
      // 跳转到目标页面，并传递数据
      const itemString = JSON.stringify(this.data);
      wx.navigateTo({
        url: `/pages/detail_edit/detail_edit?item=${itemString}`
      });
    },
    delete:function(){
      wx.navigateTo({
        url: `../testpage/testpage`
      });
      const db=wx.cloud.database({
        env:app.globalData.cloudID
      });
      db.collection("schedule").doc(this.data._id).remove({
        success: function(res) {
          wx.showToast({
            title: '删除日程成功',
            icon: 'success',
            duration: 2000
          });
        },
        fail: function(err) {
          wx.showToast({
            title: '删除失败',
            icon: 'none',
            duration: 2000
          });
        }
      })
      if(this.data.reminder==true){
      db.collection("reminder").where({
        settingValue:this.data.settingValue,
        textValue:this.data.textValue,
        titleValue:this.data.titleValue,
        timeInput:this.data.timeInput,
        urgent_level:this.data.urgent_level
      }).limit(1).get().then(res=>{
        const idToDelete = res.data[0]._id;
        db.collection('reminder').doc(idToDelete).remove();
      })
    }
  },
    deleteModal:function(){
      this.setData({
        showModalViewCancel:true
      });
    },
    cancel:function(){
      this.setData({
        showModalViewCancel:false
      });
    }
  }
})