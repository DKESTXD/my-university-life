// components/detail_edit_page/detail_edit_page.js
const app=getApp();
const date = new Date();
const years = []
const months = []
const days = []
const bigMonth = [1, 3, 5, 7, 8, 10, 12]
var getYear = date.getFullYear()
var getMonth = date.getMonth()
var getDate = date.getDate()
for (let i = getYear; i <= getYear + 40; i++) {years.push(i);}
for (let i = 1; i <= 12; i++) {months.push(i);}
for (let i = 1; i <= 31; i++) {days.push(i);}

Component({

  /**
   * 组件的属性列表
   */
  properties: {
    values:{
      type:Object,
      value:{}
    }
  },

  /**
   * 组件的初始数据
   */
  data: {
    showModalView:false,
    showModalViewCancel:false,
    hours: ["00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23"],
    minutes: ["00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59",],
    years: years,
    months: months,
    days: days,
    propDate: false,
    titleValue:"新年快乐",
    textValue:"",
    inputLength:0,
    star1:true,star2:true,star3:true,star4:true,star5:true,
    urgent_statue:"🔥🔥🔥十万火急🔥🔥🔥",
    urgent_level:"5",
    reminder:null,
    settingValue:"仅仅本次",
    timeInput:"2025-1-30-00-01",
    value:[20,0,29,"00","01"],
    hour:"00",
    minute:"01",
    year:2025,
    month:1,
    day:30
  },
  
  /**
   * 组件的方法列表
   */
  methods: {
    
    testto:function(){
      var year=this.properties.values.year;
    var month=this.properties.values.month;
    var day=this.properties.values.day;
    var hour=this.properties.values.hour;
    var minute=this.properties.values.minute;
    var timeset=`${year}-${month}-${day}T${hour}:${minute}:00`
    var timenow=`${date.getFullYear()}-${date.getMonth() + 1}-${date.getDate()}T${date.getHours()<10?`0${date.getHours()}`:date.getHours()}:${date.getMinutes()<10?`0${date.getMinutes()}`:date.getMinutes()}:${date.getSeconds()<10?`0${date.getSeconds()}`:date.getSeconds()}`

    if(timeset<=timenow&&this.data.timeInput!=''){
      wx.showModal({
        title: '截止时间似乎已经过去了哦',
        icon:'warn',
        showCancel:false
      })
    }else if(this.data.reminder==true&&this.data.timeInput==''){
      wx.showModal({
        title: '您设置了提醒，请设置截止时间哦',
        icon:'warn',
        showCancel:false
      })
    }else if(this.properties.titleValue==""){
      wx.showModal({
        title: '标题不能为空',
        icon:'warn',
        showCancel:false
      })
    }
    else{
      const db=wx.cloud.database({
        env:app.globalData.cloudID
      });
      
        if(this.properties.values.reminderold==false&&this.properties.values.reminder==true){
          db.collection("reminder").add({
            data:{
              titleValue:this.properties.values.titleValue,
              textValue:this.properties.values.textValue,
              urgent_level:this.properties.values.urgent_level,
              settingValue:this.properties.values.settingValue,
              timeInput:this.properties.values.timeInput,
              done:false
            }
          })
        }
        else if(this.properties.values.reminderold==true&&this.properties.values.reminder==false){
          db.collection("reminder").where({
            settingValue:this.properties.values.settingValueold,
            textValue:this.properties.values.textValueold,
            titleValue:this.properties.values.titleValueold,
            timeInput:this.properties.values.timeInputold,
            urgent_level:this.properties.values.urgent_levelold
          }).limit(1).get().then(res=>{
            const idToDelete = res.data[0]._id;
            db.collection('reminder').doc(idToDelete).remove();
          })
        }

      wx.navigateTo({
        url: `../testpage/testpage`
      })
      db.collection("schedule").doc(this.properties.values._id).set({
        data:{
        titleValue:this.properties.values.titleValue,
        textValue:this.properties.values.textValue,
        urgent_level:this.properties.values.urgent_level,
        reminder:this.properties.values.reminder,
        settingValue:this.properties.values.settingValue,
        timeInput:this.properties.values.timeInput,
        setTime:this.properties.values.setTime,
        },
        success:res=>{
          wx.showToast({
            title: '日程编辑成功',
          })
          console.log("日程编辑成功",res._id)
        },
        fail:err=>{
          wx.showToast({
            icon:noneParamsEaseFuncs,
            title: '日程编辑失败',
          })
          console.log("日程编辑失败",err)
        }
      });


    }
    },
    titleInput: function(e) {
      const newTitle = e.detail.value;
      this.triggerEvent('updatevalues', { key: 'titleValue', value: newTitle });
    },
    textInput: function(e) {
      const newText = e.detail.value;
      const newlength= newText.length;
      this.triggerEvent('updatevalues', { key: 'textValue', value: newText });
      this.triggerEvent('updatevalues', { key: 'inputLength', value: newlength });
    },
    star1_tap: function() {
      this.setData({
        star1: true, star2: false, star3: false, star4: false, star5: false, urgent_statue: '以后的事以后再说',
        urgent_level: 1
      });

      this.triggerEvent('updatevalues', { key: 'urgent_statue', value: this.data.urgent_statue });
      this.triggerEvent('updatevalues', { key: 'urgent_level',value: this.data.urgent_level });
      this.triggerEvent('updatevalues', { key: 'star1', value: this.data.star1 });
      this.triggerEvent('updatevalues', { key: 'star2', value: this.data.star2 });
      this.triggerEvent('updatevalues', { key: 'star3', value: this.data.star3 });
      this.triggerEvent('updatevalues', { key: 'star4', value: this.data.star4 });
      this.triggerEvent('updatevalues', { key: 'star5', value: this.data.star5 });
    },
    star2_tap: function() {
      this.setData({
        star1: true, star2: true, star3: false, star4: false, star5: false, urgent_statue: '轻轻松 时间还很充足',
        urgent_level: 2
      });
      this.triggerEvent('updatevalues', { key: 'urgent_statue', value: this.data.urgent_statue });
      this.triggerEvent('updatevalues', { key: 'urgent_level',value: this.data.urgent_level });
      this.triggerEvent('updatevalues', { key: 'star1', value: this.data.star1 });
      this.triggerEvent('updatevalues', { key: 'star2', value: this.data.star2 });
      this.triggerEvent('updatevalues', { key: 'star3', value: this.data.star3 });
      this.triggerEvent('updatevalues', { key: 'star4', value: this.data.star4 });
      this.triggerEvent('updatevalues', { key: 'star5', value: this.data.star5 });
    },
    star3_tap: function() {
      this.setData({
        star1: true, star2: true, star3: true, star4: false, star5: false, urgent_statue: '你意识到你不能再拖了',
        urgent_level: 3
      });
      this.triggerEvent('updatevalues', { key: 'urgent_statue', value: this.data.urgent_statue });
      this.triggerEvent('updatevalues', { key: 'urgent_level',value: this.data.urgent_level });
      this.triggerEvent('updatevalues', { key: 'star1', value: this.data.star1 });
      this.triggerEvent('updatevalues', { key: 'star2', value: this.data.star2 });
      this.triggerEvent('updatevalues', { key: 'star3', value: this.data.star3 });
      this.triggerEvent('updatevalues', { key: 'star4', value: this.data.star4 });
      this.triggerEvent('updatevalues', { key: 'star5', value: this.data.star5 });
    },
    star4_tap: function() {
      this.setData({
        star1: true, star2: true, star3: true, star4: true, star5: false, urgent_statue: '非常重要 但是似乎会有比它更重要的事',
        urgent_level: 4
      });
      this.triggerEvent('updatevalues', { key: 'urgent_statue', value: this.data.urgent_statue });
      this.triggerEvent('updatevalues', { key: 'urgent_level',value: this.data.urgent_level });
      this.triggerEvent('updatevalues', { key: 'star1', value: this.data.star1 });
      this.triggerEvent('updatevalues', { key: 'star2', value: this.data.star2 });
      this.triggerEvent('updatevalues', { key: 'star3', value: this.data.star3 });
      this.triggerEvent('updatevalues', { key: 'star4', value: this.data.star4 });
      this.triggerEvent('updatevalues', { key: 'star5', value: this.data.star5 });
    },
    star5_tap: function() {
      this.setData({
        star1: true, star2: true, star3: true, star4: true, star5: true, urgent_statue: '🔥🔥🔥十万火急🔥🔥🔥',
        urgent_level: 5
      });
      this.triggerEvent('updatevalues', { key: 'urgent_statue', value: this.data.urgent_statue });
      this.triggerEvent('updatevalues', { key: 'urgent_level',value: this.data.urgent_level });
      this.triggerEvent('updatevalues', { key: 'star1', value: this.data.star1 });
      this.triggerEvent('updatevalues', { key: 'star2', value: this.data.star2 });
      this.triggerEvent('updatevalues', { key: 'star3', value: this.data.star3 });
      this.triggerEvent('updatevalues', { key: 'star4', value: this.data.star4 });
      this.triggerEvent('updatevalues', { key: 'star5', value: this.data.star5 });
    },
    reminder_change: function(e) {
      const value1 = e.detail.value; 
    const that = this; 
    if (value1 === true) {
      this.setData({
        reminder: true
    });
    this.triggerEvent('updatevalues', { key: 'reminder', value: this.data.reminder });
        wx.requestSubscribeMessage({
            tmplIds: ['v3J9x3kz6FYW3Xt4XwJ-d3qpbTU-J19XBh91HwOKqLA'],
            success(res) {
                if (res['v3J9x3kz6FYW3Xt4XwJ-d3qpbTU-J19XBh91HwOKqLA'] === 'accept') {
                    wx.showToast({
                        title: '成功订阅日程提醒',
                        icon: "none"
                    });
                } else if (res['v3J9x3kz6FYW3Xt4XwJ-d3qpbTU-J19XBh91HwOKqLA'] === 'reject') {
                    wx.showToast({
                        title: '您将无法使用日程提醒服务',
                        icon: "none"
                    });
                    that.setData({
                        reminder: false 
                    });
                    that.triggerEvent('updatevalues', { key: 'reminder', value: that.data.reminder });
                }
            },
            fail(err) {
                console.error('订阅消息调用失败', err);
                wx.showToast({
                    title: '订阅失败，请检查设置',
                    icon: "none"
                });
                that.setData({
                    reminder: false 
                });
                that.triggerEvent('updatevalues', { key: 'reminder', value: that.data.reminder });
            }
        });
    } else {
        this.setData({
            reminder: false
        });
        this.triggerEvent('updatevalues', { key: 'reminder', value: this.data.reminder });
    }
    },
    showModal: function() {
      this.setData({
        showModalView: true
      });
    },
    selectSetting_only: function(e) {
      this.setData({
        settingValue: "仅仅本次",
        showModalView: false
      });
      this.triggerEvent('updatevalues', { key: 'settingValue', value: this.data.settingValue });
    },
    selectSetting_weekly: function(e) {
      this.setData({
        settingValue: "每周一次",
        showModalView: false
      });
      this.triggerEvent('updatevalues', { key: 'settingValue', value: this.data.settingValue });
    },
    selectSetting_monthly: function(e) {
      this.setData({
        settingValue: "每月一次",
        showModalView: false
      });
      this.triggerEvent('updatevalues', { key: 'settingValue', value: this.data.settingValue });
    },
    cancel: function() {
      this.setData({
        showModalViewCancel: false
      });
    },
    cancel_tap: function() {
      this.setData({
        showModalViewCancel: true
      });
    },
    confirm_edit: function() {
      this.setData({
        showModalViewCancel: false
      });
      wx.navigateTo({
        url: `../testpage/testpage`
      });
    },
    closePick: function() {
      this.setData({
        propDate: false
      });
    },
    openPick: function() {
      let valueTime = '';
      if (this.data.values.timeInput) {
        const res = this.data.values.timeInput.split("-");
        const year = res[0].substring(2);
        const newRes = [year - 25, res[1] - 1, res[2] - 1, res[3], res[4]];
        valueTime = newRes;
      } else {
        valueTime = [this.data.year - 25, this.data.month - 1, this.data.day - 1, "00", "00"];
      }
      this.setData({
        propDate: true,
        value: valueTime,
      });
      this.triggerEvent('updatevalues', { key: 'value', value: this.data.value });
    },
    clearPick: function() {
      this.setData({
        timeInput: '',
        value: '',
      });
      this.triggerEvent('updatevalues', { key: 'timeInput', value: this.data.timeInput });
      this.triggerEvent('updatevalues', { key: 'value', value: this.data.value });
    },
    contains: function(arr, obj) {
      var i = arr.length;
      while (i--) {
        if (arr[i] === obj) {
          return true;
        }
      }
      return false;
    },
    setDays: function(day) {
      const temp = [];
      for (let i = 1; i <= day; i++) {
        temp.push(i);
      }
      this.setData({
        days: temp,
      });
    },
    bindChange: function(e) {
      const val = e.detail.value;
      const setYear = this.data.years[val[0]];
      const setMonth = this.data.months[val[1]];
      const setDay = this.data.days[val[2]];
      const setHour = this.data.hours[val[3]];
      const setMinute = this.data.minutes[val[4]];
      this.setData({
        year: setYear,
        month: setMonth,
        day: setDay,
        hour: setHour,
        minute: setMinute,
      });
      if (setMonth === 2) {
        this.setDays((setYear % 4 === 0 && setYear % 100 !== 0) ? 29 : 28);
      } else {
        this.setDays(this.contains(bigMonth, setMonth) ? 31 : 30);
      }
      this.setData({
        timeInput: setYear + '-' + setMonth + '-' + setDay + '-' + setHour + '-' + setMinute
      });
      this.triggerEvent('updatevalues', { key: 'timeInput', value: this.data.timeInput });
      this.triggerEvent('updatevalues', { key: 'year', value: this.data.year });
      this.triggerEvent('updatevalues', { key: 'month', value: this.data.month });
      this.triggerEvent('updatevalues', { key: 'day', value: this.data.day });
      this.triggerEvent('updatevalues', { key: 'hour', value: this.data.hour });
      this.triggerEvent('updatevalues', { key: 'minute', value: this.data.minute });
      this.openPick();
    },
  }
})